"""Part I — fibre channel model and distance-driven SNS-TF-QKD key rate.

Corresponds to notebook Sections 3 (Channel Model) and 9 (Parameter
Optimisation), plus the fibre-specific wrapper used throughout Part I
(Sections 8-12) around the canonical eta-native engine in ``src/engine.py``.
"""

import numpy as np
from scipy.optimize import minimize

from src.engine import (
    ETA_D, P_DARK, F_EC,
    gain_X, error_X, S_Z_gain, E_Z_qber,
    s1_decoy_lower_bound, e1ph_upper_bound, binary_entropy,
)

# ---------------------------------------------------------------------------
# Section 3: Channel model
# ---------------------------------------------------------------------------
ALPHA = 0.2      # dB/km, standard telecom fibre attenuation


def channel_eta(L_km, alpha=ALPHA, eta_d=ETA_D):
    """Single-photon transmittance from Alice (or Bob) to one of Charlie's
    detectors, for a total Alice-Bob distance L_km, assuming Charlie sits at
    the mid-point (each arm has length L_km/2).

        eta(L) = eta_d * 10^(-alpha * (L/2) / 10)

    This is the key quantity responsible for TF-QKD's square-root distance
    scaling: eta(L) ~ 10^(-alpha L / 20), i.e. the square root of the
    transmittance 10^(-alpha L /10) of the *full* distance L that governs
    ordinary point-to-point QKD.
    """
    d_arm = L_km / 2.0
    eta_channel = 10.0 ** (-alpha * d_arm / 10.0)
    return eta_d * eta_channel


def true_single_photon_yield(eta, p_dark=P_DARK):
    """Yield of an X1-window: click probability given exactly one (total)
    photon is incident from either side onto the beamsplitter with single-arm
    transmittance eta. Used only for the Section 6 validation."""
    return 1.0 - (1.0 - 2.0 * p_dark) * (1.0 - eta)


# ---------------------------------------------------------------------------
# Section 8: Eq. (4) - asymptotic key rate per time window (fibre wrapper)
# ---------------------------------------------------------------------------
def sns_key_rate(eps, mu_prime, L, e_a, decoy_ratio=0.2, f_ec=F_EC, p_dark=P_DARK):
    """Full asymptotic SNS-TFQKD key rate per time window, Eq. (4), for
    Alice-Bob distance L (km) and single-photon-interference misalignment
    error e_a.

    Decoy intensities: mu2 = mu_prime (matched to the signal intensity),
                        mu1 = decoy_ratio * mu2 (a weaker decoy), mu0 = 0
                        (vacuum).
    This is a standard, documented choice; the paper does not fix specific
    decoy intensities for the asymptotic (infinite-decoy) simulation.

    Returns the key rate R (bits per pulse); may be negative if no secure key
    can be extracted with the given parameters.
    """
    eta = channel_eta(L)
    mu2 = mu_prime
    mu1 = decoy_ratio * mu2
    s0 = 2.0 * p_dark

    S1 = gain_X(mu1, eta, p_dark)
    S2 = gain_X(mu2, eta, p_dark)
    E1 = error_X(mu1, eta, e_a, p_dark)

    s1 = s1_decoy_lower_bound(mu1, mu2, S1, S2, s0)
    s1 = max(s1, 1e-15)

    e1ph = e1ph_upper_bound(mu1, S1, E1, s0, s1)
    e1ph = float(np.clip(e1ph, 0.0, 0.5))

    Sz = S_Z_gain(eps, mu_prime, eta, p_dark)
    Ez = E_Z_qber(eps, mu_prime, eta, p_dark)
    Ez = float(np.clip(Ez, 1e-12, 0.5))

    R = (2 * eps * (1 - eps) * mu_prime * np.exp(-mu_prime) * s1 * (1 - binary_entropy(e1ph))
         - Sz * f_ec * binary_entropy(Ez))
    return R


# ---------------------------------------------------------------------------
# Section 9: numerical optimisation of (eps, mu') for maximum key rate
# ---------------------------------------------------------------------------
def optimize_key_rate(L, e_a, decoy_ratio=0.2,
                       eps_bounds=(1e-5, 0.5), mu_bounds=(1e-5, 1.0),
                       n_grid=25):
    """Two-stage optimisation (coarse grid + Nelder-Mead refinement) of the
    SNS-TFQKD key rate Eq. (4) over the sending probability eps and signal
    intensity mu'.

    Returns (R_opt, eps_opt, mu_opt).
    """
    eps_grid = np.logspace(np.log10(eps_bounds[0]), np.log10(eps_bounds[1]), n_grid)
    mu_grid = np.logspace(np.log10(mu_bounds[0]), np.log10(mu_bounds[1]), n_grid)

    best_R = -np.inf
    best_pt = (eps_grid[0], mu_grid[0])
    for e in eps_grid:
        for m in mu_grid:
            R = sns_key_rate(e, m, L, e_a, decoy_ratio)
            if R > best_R:
                best_R = R
                best_pt = (e, m)

    def negative_rate(x):
        e, m = x
        if not (eps_bounds[0] < e < eps_bounds[1]) or not (mu_bounds[0] < m < mu_bounds[1]):
            return 1e3
        return -sns_key_rate(e, m, L, e_a, decoy_ratio)

    res = minimize(negative_rate, x0=list(best_pt), method='Nelder-Mead',
                    options={'xatol': 1e-9, 'fatol': 1e-13, 'maxiter': 800, 'maxfev': 800})

    if res.success and -res.fun > best_R:
        return -res.fun, res.x[0], res.x[1]
    return best_R, best_pt[0], best_pt[1]
