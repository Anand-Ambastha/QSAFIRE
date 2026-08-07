"""Part 1 — Channel Models (deterministic FSO path loss + Gamma-Gamma
atmospheric turbulence), and Part 2 — turbulence-regime classification.
"""

import numpy as np
from scipy.special import gamma as gamma_fn, kv as bessel_k

from src.fso.config import ALPHA_FSO_DB_KM, K_WAVENUMBER

# ---------------------------------------------------------------------------
# Model A: deterministic FSO path loss
# ---------------------------------------------------------------------------
def fso_deterministic_eta(d_km, alpha=ALPHA_FSO_DB_KM, eta_d=0.8):
    """Model A - deterministic FSO channel transmissivity (single arm):
        eta(d) = eta_d * 10^(-alpha*d/10)
    d_km: one-way propagation distance (km); alpha: dB/km; eta_d: detector
    efficiency.
    """
    return eta_d * 10.0 ** (-alpha * np.asarray(d_km, dtype=float) / 10.0)


def channel_loss_dB(d_km, alpha=ALPHA_FSO_DB_KM):
    """Pure atmospheric channel loss in dB (excludes detector efficiency)."""
    return alpha * np.asarray(d_km, dtype=float)


# ---------------------------------------------------------------------------
# Model B: Gamma-Gamma atmospheric turbulence
# ---------------------------------------------------------------------------
def rytov_variance(d_km, Cn2, k=K_WAVENUMBER):
    """Plane-wave Rytov variance:  sigma_R^2 = 1.23 * Cn2 * k^(7/6) * L^(11/6)
    d_km: one-way propagation distance (km) -> converted to metres
    internally.
    """
    L_m = np.asarray(d_km, dtype=float) * 1000.0
    return 1.23 * Cn2 * k ** (7.0 / 6.0) * L_m ** (11.0 / 6.0)


SIGMA_R2_VALIDITY_MAX = 25.0
"""Ceiling on the plane-wave Rytov variance for which the moment-based
alpha_g/beta_g formula below is trustworthy.

BUG NOTE (found and fixed): for sigma_R^2 well beyond this range, the
denominators 1.11*sigma_R^(12/5) and 0.69*sigma_R^(12/5) grow faster than
the sigma_R^2 numerator (exponents 12/5 * 7/6 = 2.8 and 12/5 * 5/6 = 2.0,
both > 1), so the exponent argument -> 0 and BOTH alpha_g and beta_g diverge
to very large values instead of the physically correct saturation-regime
behaviour (alpha_g -> infinity while beta_g -> ~1, keeping the scintillation
index sigma_I^2 = 1/alpha_g + 1/beta_g + 1/(alpha_g*beta_g) near 1).
Concretely: sigma_I^2 peaks near sigma_R^2~1 (~0.71) and then falls back
toward 0 as sigma_R^2 grows into the hundreds/thousands - backwards from the
known saturation phenomenon. Once alpha_g/beta_g are that large,
sample_gamma_gamma() draws h with ~zero variance (h~1 almost surely), so the
Gamma-Gamma channel eta_i = eta_det*h collapses numerically onto the
deterministic channel eta_det. This is why 'gg' and 'deterministic' results
were coming out identical: at the distances used (DISTANCES_MAIN reaches
sigma_R^2 in the thousands by 200 km), the model was being extrapolated far
outside where the formula is valid (weak/moderate turbulence, roughly
sigma_R^2 <~ 10-25 - see classify_turbulence_regime's own 'strong' caveat,
which flagged this but didn't act on it). Fixing this properly (a strong-
fluctuation/saturation-regime parameterisation) is out of scope here; instead
downstream code must not treat gamma_gamma_params()/sample_gamma_gamma()
output as physically meaningful above SIGMA_R2_VALIDITY_MAX - see
sigma_R2_within_validity() and its use in metrics.build_gamma_gamma_table.
"""


def sigma_R2_within_validity(sigma_R2, ceiling=SIGMA_R2_VALIDITY_MAX):
    """True where the weak/moderate-turbulence alpha_g/beta_g formula is
    still within its documented validity range (see SIGMA_R2_VALIDITY_MAX)."""
    return np.asarray(sigma_R2, dtype=float) <= ceiling


def gamma_gamma_params(sigma_R2):
    """Large-scale (alpha_g) and small-scale (beta_g) Gamma-Gamma shape
    parameters from the plane-wave Rytov variance:

        alpha_g = [exp(0.49*sigma_R^2 / (1+1.11*sigma_R^(12/5))^(7/6)) - 1]^-1
        beta_g  = [exp(0.51*sigma_R^2 / (1+0.69*sigma_R^(12/5))^(5/6)) - 1]^-1

    NOTE: only trustworthy for sigma_R^2 <~ SIGMA_R2_VALIDITY_MAX (see that
    constant's docstring) - beyond this the formula is being extrapolated out
    of its valid domain and alpha_g/beta_g diverge unphysically. Callers that
    aggregate over a distance range spanning both regimes (e.g.
    metrics.build_gamma_gamma_table) must filter on
    sigma_R2_within_validity() before drawing conclusions from the output.
    """
    s2 = np.asarray(sigma_R2, dtype=float)
    alpha_g = 1.0 / (np.expm1(0.49 * s2 / (1.0 + 1.11 * s2 ** (12.0 / 5.0)) ** (7.0 / 6.0)))
    beta_g = 1.0 / (np.expm1(0.51 * s2 / (1.0 + 0.69 * s2 ** (12.0 / 5.0)) ** (5.0 / 6.0)))
    return alpha_g, beta_g


def classify_turbulence_regime(sigma_R2):
    """Part 2 validity check. Returns one of:
      'weak'      : sigma_R^2 < 1           -> Gamma-Gamma valid (reduces
                                                 towards log-normal)
      'moderate'  : 1 <= sigma_R^2 < 10      -> Gamma-Gamma valid (its main
                                                 regime of superiority)
      'strong'    : sigma_R^2 >= 10          -> approaching/within saturation;
                                                 Gamma-Gamma still used but
                                                 flagged, since irradiance
                                                 moments derived under the
                                                 Rytov approximation become
                                                 less reliable without
                                                 aperture-averaging
                                                 corrections.
    """
    s2 = np.asarray(sigma_R2, dtype=float)
    regime = np.where(s2 < 1.0, "weak", np.where(s2 < 10.0, "moderate", "strong"))
    return regime


def validity_flag(sigma_R2):
    """Human-readable validity flag per Part 2:
      sigma_R^2 < 1   -> 'Gamma-Gamma valid (weak-fluctuation regime;
                          log-normal also acceptable)'
      1<=sigma_R^2<10  -> 'Gamma-Gamma valid (moderate turbulence, GG''s
                          primary regime)'
      sigma_R^2 >= 10  -> 'Saturation regime - Gamma-Gamma applied with
                          caution'
    """
    s2 = np.asarray(sigma_R2, dtype=float)
    out = []
    it = np.atleast_1d(s2)
    for v in it:
        if v < 1.0:
            out.append("Gamma-Gamma valid / log-normal preferred (weak fluctuations)")
        elif v < 10.0:
            out.append("Gamma-Gamma valid (moderate turbulence)")
        else:
            out.append("Saturation regime - Gamma-Gamma applied with caution")
    result = np.array(out, dtype=object)
    return result if s2.ndim else result[0]


def gamma_gamma_pdf(h, alpha_g, beta_g):
    """Gamma-Gamma PDF (for validation / plotting only - sampling itself uses
    the Gamma-product method below, which is exact and does not rely on
    evaluating this PDF):

        f(h) = 2*(alpha_g*beta_g)^((alpha_g+beta_g)/2) / (Gamma(alpha_g)*Gamma(beta_g))
               * h^((alpha_g+beta_g)/2 - 1) * K_{alpha_g-beta_g}(2*sqrt(alpha_g*beta_g*h))
    """
    h = np.asarray(h, dtype=float)
    pref = 2.0 * (alpha_g * beta_g) ** ((alpha_g + beta_g) / 2.0) / (gamma_fn(alpha_g) * gamma_fn(beta_g))
    with np.errstate(divide="ignore", invalid="ignore"):
        val = pref * h ** ((alpha_g + beta_g) / 2.0 - 1.0) * bessel_k(alpha_g - beta_g, 2.0 * np.sqrt(alpha_g * beta_g * h))
    return np.nan_to_num(val, nan=0.0, posinf=0.0, neginf=0.0)


def sample_gamma_gamma(alpha_g, beta_g, size, rng):
    """Exact Gamma-Gamma sampler via the standard doubly-stochastic
    construction:
    h = X*Y, X~Gamma(shape=alpha_g, scale=1/alpha_g) [unit mean],
             Y~Gamma(shape=beta_g,  scale=1/beta_g)  [unit mean]
    so that E[h] = E[X]*E[Y] = 1 (unit-mean fading, as required so that the
    ENSEMBLE mean channel equals the deterministic Model-A channel).
    """
    X = rng.gamma(shape=alpha_g, scale=1.0 / alpha_g, size=size)
    Y = rng.gamma(shape=beta_g, scale=1.0 / beta_g, size=size)
    return X * Y