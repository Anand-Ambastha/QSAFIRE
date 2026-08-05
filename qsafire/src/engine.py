"""SNS-TF-QKD physics engine (Wang, Yu & Hu, arXiv:1805.09222 / PRA 98, 062323 (2018)).

This module is the **single canonical source** of the protocol's core equations
(photon-number weights Eqs. 6/8/9, gain/error model Section 4, Z-window gain/QBER
Section 5.3, decoy-state single-photon-yield bound Eq. (44), phase-error upper
bound Eq. (45), binary entropy, and the asymptotic key-rate formula Eq. (4)),
expressed directly in terms of the channel transmittance ``eta``.

Origin note: the source notebook defines this exact set of functions **twice**,
verbatim, in two different narrative sections:

- Sections 5-8 (Part I, distance-driven, fibre channel): the functions operate
  on a fibre distance ``L`` via ``channel_eta(L)`` (see
  ``src/simulation/fiber_channel.py``).
- Section 14 / "Part II" (FSO/turbulence extension): the notebook explicitly
  states this is *"a verbatim re-statement of Sections 5-8 above ... with the
  single architectural change"* that every function accepts ``eta`` directly,
  so the same equations can be driven by the deterministic FSO path-loss model,
  Monte Carlo turbulence samples, or Gauss-Laguerre quadrature nodes.

Because both notebook cells are byte-for-byte identical apart from that one
signature detail, they are consolidated here into one module rather than
duplicated, exactly matching the notebook's own description of the second copy
as a "verbatim restatement". Both ``src/simulation/fiber_channel.py`` (Part I)
and ``src/fso`` (Part II) import from here, so behaviour is unchanged and no
computation is duplicated or lost.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Fixed detector / protocol parameters (Section III of Wang, Yu & Hu 2018)
# ---------------------------------------------------------------------------
ETA_D = 0.8        # detector efficiency
P_DARK = 1e-11     # dark-count probability per detector per pulse
F_EC = 1.1         # error-correction inefficiency factor
DECOY_RATIO = 0.2  # mu1 = DECOY_RATIO * mu2  (weak decoy relative to signal-matched decoy)


def p0_weight(mu):
    """Vacuum photon-number weight p0(mu) = exp(-2*mu). Eq. (6)."""
    return np.exp(-2.0 * mu)


def p1_weight(mu):
    """Single-photon weight p1(mu) = 2*mu*exp(-2*mu). Eq. (8)."""
    return 2.0 * mu * np.exp(-2.0 * mu)


def p2_weight(mu):
    """Two-photon weight p2(mu) = 2*mu^2*exp(-2*mu). Eq. (9)."""
    return 2.0 * mu ** 2 * np.exp(-2.0 * mu)


def s0_vacuum_yield(p_dark=P_DARK):
    """Vacuum yield s0: dark-count-only effective-event probability."""
    return 2.0 * p_dark


def gain_X(mu, eta, p_dark=P_DARK):
    """Observed click rate (gain) S_mu for an X-window (decoy window) in which
    BOTH Alice and Bob send a coherent state of the SAME intensity mu.

    Total mean photon number arriving at Charlie (summed over both arms) is
    2*mu*eta. S_mu = 1 - (1 - 2*p_dark) * exp(-2*mu*eta)
    """
    return 1.0 - (1.0 - 2.0 * p_dark) * np.exp(-2.0 * mu * eta)


def error_X(mu, eta, e_a, p_dark=P_DARK):
    """X-basis (interference) error rate E_mu^X for decoy intensity mu.

    - The "correct-photon" contribution 1 - exp(-2*mu*eta) suffers the
      misalignment error rate e_a.
    - The dark-count contribution 2*p_dark*exp(-2*mu*eta) is unbiased noise
      -> error 1/2.
    """
    correct_click_prob = 1.0 - np.exp(-2.0 * mu * eta)
    dark_click_prob = 2.0 * p_dark * np.exp(-2.0 * mu * eta)
    S = gain_X(mu, eta, p_dark)
    numerator = e_a * correct_click_prob + 0.5 * dark_click_prob
    return np.where(S > 0, numerator / np.where(S > 0, S, 1.0), 0.5)


def S_Z_gain(eps, mu_prime, eta, p_dark=P_DARK):
    """Observed Z-window effective-event rate (gain). See Section 5.3."""
    one_side_click = 1.0 - np.exp(-mu_prime * eta)
    both_send_click = 1.0 - np.exp(-2.0 * mu_prime * eta)
    neither_send_click = 2.0 * p_dark
    return (2.0 * eps * (1.0 - eps) * one_side_click
            + eps ** 2 * both_send_click
            + (1.0 - eps) ** 2 * neither_send_click)


def E_Z_qber(eps, mu_prime, eta, p_dark=P_DARK):
    """Z-basis quantum bit error rate. See Section 5.3."""
    both_send_click = 1.0 - np.exp(-2.0 * mu_prime * eta)
    neither_send_click = 2.0 * p_dark
    Sz = S_Z_gain(eps, mu_prime, eta, p_dark)
    numerator = eps ** 2 * both_send_click + (1.0 - eps) ** 2 * neither_send_click
    return np.where(Sz > 0, numerator / np.where(Sz > 0, Sz, 1.0), 0.5)


def s1_decoy_lower_bound(mu1, mu2, S_mu1, S_mu2, s0):
    """Equation (44): 3-intensity decoy-state lower bound on the single-photon
    yield s1, given mu0 = 0."""
    denominator = p2_weight(mu2) * p1_weight(mu1) - p2_weight(mu1) * p1_weight(mu2)
    numerator = (p2_weight(mu2) * (S_mu1 - p0_weight(mu1) * s0)
                 - p2_weight(mu1) * (S_mu2 - p0_weight(mu2) * s0))
    return numerator / denominator


def e1ph_upper_bound(mu1, S_mu1, E_mu1, s0, s1):
    """Equation (45): upper bound on the single-photon phase-flip error rate
    e1^ph, obtained from the X-basis observed error rate at intensity mu1."""
    numerator = S_mu1 * E_mu1 - np.exp(-2.0 * mu1) * s0 / 2.0
    denominator = 2.0 * mu1 * np.exp(-2.0 * mu1) * s1
    return numerator / denominator


def binary_entropy(x):
    """H(x) = -x log2 x - (1-x) log2 (1-x), safely clipped to (0,1)."""
    x = np.clip(x, 1e-15, 1 - 1e-15)
    return -x * np.log2(x) - (1 - x) * np.log2(1 - x)


def sns_key_rate_eta(eps, mu_prime, eta, e_a, decoy_ratio=DECOY_RATIO,
                      f_ec=F_EC, p_dark=P_DARK):
    """Eq. (4) of Wang, Yu & Hu (2018), evaluated at an arbitrary, directly
    supplied channel transmittance `eta` (single-arm, Alice/Bob -> Charlie).

    Fully vectorised: eps, mu_prime, eta, e_a may be scalars or numpy arrays
    that broadcast together (this is what makes the Monte Carlo module fast:
    N_MC = 1e5 turbulence realizations are evaluated in one call, no Python
    loop).

    Returns key rate R (bits/pulse), which may be negative (no secure key).
    """
    eta = np.asarray(eta, dtype=float)
    mu2 = mu_prime
    mu1 = decoy_ratio * mu2
    s0 = s0_vacuum_yield(p_dark)

    S1 = gain_X(mu1, eta, p_dark)
    S2 = gain_X(mu2, eta, p_dark)
    E1 = error_X(mu1, eta, e_a, p_dark)

    s1 = s1_decoy_lower_bound(mu1, mu2, S1, S2, s0)
    s1 = np.maximum(s1, 1e-15)

    e1ph = e1ph_upper_bound(mu1, S1, E1, s0, s1)
    e1ph = np.clip(e1ph, 0.0, 0.5)

    Sz = S_Z_gain(eps, mu_prime, eta, p_dark)
    Ez = E_Z_qber(eps, mu_prime, eta, p_dark)
    Ez = np.clip(Ez, 1e-12, 0.5)

    R = (2 * eps * (1 - eps) * mu_prime * np.exp(-mu_prime) * s1 * (1 - binary_entropy(e1ph))
         - Sz * f_ec * binary_entropy(Ez))
    return R


def full_metrics_eta_vec(eps, mu_prime, eta, e_a, decoy_ratio=DECOY_RATIO,
                          f_ec=F_EC, p_dark=P_DARK):
    """Vectorised sibling of full_metrics_eta: eta may be an array (e.g.
    Gauss-Laguerre nodes or Monte Carlo samples). Returns a dict of arrays,
    same length as eta."""
    eta = np.asarray(eta, dtype=float)
    mu2 = mu_prime
    mu1 = decoy_ratio * mu2
    s0 = s0_vacuum_yield(p_dark)
    S1 = gain_X(mu1, eta, p_dark)
    S2 = gain_X(mu2, eta, p_dark)
    E1 = error_X(mu1, eta, e_a, p_dark)
    s1 = np.maximum(s1_decoy_lower_bound(mu1, mu2, S1, S2, s0), 1e-15)
    e1ph = np.clip(e1ph_upper_bound(mu1, S1, E1, s0, s1), 0.0, 0.5)
    Sz = S_Z_gain(eps, mu_prime, eta, p_dark)
    Ez = np.clip(E_Z_qber(eps, mu_prime, eta, p_dark), 1e-12, 0.5)
    R = sns_key_rate_eta(eps, mu_prime, eta, e_a, decoy_ratio, f_ec, p_dark)
    return dict(eta=eta, Qxx=S2, Exx=E1, s1=s1, e1ph=e1ph, Sz=Sz, Ez=Ez, R=R)


def full_metrics_eta(eps, mu_prime, eta, e_a, decoy_ratio=DECOY_RATIO,
                      f_ec=F_EC, p_dark=P_DARK):
    """Return a dict of every intermediate quantity (Part 8 performance
    metrics), for a given eta. Scalar-oriented (used for the deterministic-
    channel table)."""
    eta = float(eta)
    mu2 = mu_prime
    mu1 = decoy_ratio * mu2
    s0 = s0_vacuum_yield(p_dark)
    S1 = gain_X(mu1, eta, p_dark)
    S2 = gain_X(mu2, eta, p_dark)
    E1 = error_X(mu1, eta, e_a, p_dark)
    s1 = max(float(s1_decoy_lower_bound(mu1, mu2, S1, S2, s0)), 1e-15)
    e1ph = float(np.clip(e1ph_upper_bound(mu1, S1, E1, s0, s1), 0.0, 0.5))
    Sz = float(S_Z_gain(eps, mu_prime, eta, p_dark))
    Ez = float(np.clip(E_Z_qber(eps, mu_prime, eta, p_dark), 1e-12, 0.5))
    R = float(sns_key_rate_eta(eps, mu_prime, eta, e_a, decoy_ratio, f_ec, p_dark))
    return dict(eta=eta, Qxx=S2, Exx=E1, s1=s1, e1ph=e1ph, Sz=Sz, Ez=Ez, R=R)
