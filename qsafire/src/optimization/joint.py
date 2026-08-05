"""Joint signal/decoy intensity and sending-probability optimization
(notebook Part X / Phase 8).

Fixes two issues left in place through Phases 1-7: (1) epsilon was bounded to
(0.01, 0.99) in every downstream optimizer wrapper, when Part I's own
documented convention is epsilon in (0, 0.5) (send/don't-send label
symmetry); (2) the decoy intensity mu1 was always fixed at a flat
DECOY_RATIO=0.2*mu2 rather than independently optimized. Both are corrected
here via a new rate function with mu1 as a free argument, and joint 3-D
optimizers (grid + Nelder-Mead) replacing the Part V/Part VIII wrappers.

Origin: notebook cells 198 (sns_key_rate_eta_joint + refactor-correctness
check), 200 (optimize_joint_signal_decoy), 202 (optimize_joint_ensemble).
"""

import numpy as np
from scipy.optimize import minimize

from src.engine import (
    F_EC, P_DARK, s0_vacuum_yield, gain_X, error_X, s1_decoy_lower_bound,
    e1ph_upper_bound, S_Z_gain, E_Z_qber, binary_entropy, sns_key_rate_eta,
)
from src.fso.channel_models import gamma_gamma_params, sample_gamma_gamma


def sns_key_rate_eta_joint(eps, mu2, mu1, eta, e_a, f_ec=F_EC, p_dark=P_DARK):
    """Identical physics to sns_key_rate_eta (Part I), but mu1 (decoy
    intensity) is now a free argument instead of being derived from a fixed
    decoy_ratio. Cell 198."""
    eta = np.asarray(eta, dtype=float)
    s0 = s0_vacuum_yield(p_dark)
    S1 = gain_X(mu1, eta, p_dark)
    S2 = gain_X(mu2, eta, p_dark)
    E1 = error_X(mu1, eta, e_a, p_dark)
    s1 = np.maximum(s1_decoy_lower_bound(mu1, mu2, S1, S2, s0), 1e-15)
    e1ph = np.clip(e1ph_upper_bound(mu1, S1, E1, s0, s1), 0.0, 0.5)
    Sz = S_Z_gain(eps, mu2, eta, p_dark)
    Ez = np.clip(E_Z_qber(eps, mu2, eta, p_dark), 1e-12, 0.5)
    R = (2*eps*(1-eps)*mu2*np.exp(-mu2)*s1*(1-binary_entropy(e1ph)) - Sz*f_ec*binary_entropy(Ez))
    return R


def refactor_correctness_check():
    """Cell 198 (script part): fixing mu1 = 0.2*mu2 must exactly reproduce
    the original sns_key_rate_eta."""
    _eta_chk = 0.01
    _max_diff = 0.0
    for _eps_t, _mu_t in [(0.02, 0.2), (0.05, 0.5), (0.1, 0.1), (0.3, 0.8)]:
        _R_old = float(sns_key_rate_eta(_eps_t, _mu_t, _eta_chk, 0.15, decoy_ratio=0.2))
        _R_new = float(sns_key_rate_eta_joint(_eps_t, _mu_t, 0.2*_mu_t, _eta_chk, 0.15))
        _max_diff = max(_max_diff, abs(_R_old - _R_new))
    print(f"Max difference vs. original sns_key_rate_eta at decoy_ratio=0.2 (refactor check): {_max_diff:.2e}")
    assert _max_diff < 1e-12, "Refactored joint rate function does not match the original at fixed decoy ratio"
    print("PASS: refactor is exact.")
    return _max_diff


def optimize_joint_signal_decoy(eta_value, e_a, p_dark, n_grid=14,
                                 eps_bounds=(1e-5, 0.5), mu2_bounds=(1e-3, 1.0), frac_bounds=(0.005, 0.95)):
    """Log-spaced coarse grid (Part I convention) + Nelder-Mead refinement.
    mu1 = frac*mu2 with frac in (0,1) so mu1 < mu2 is automatically enforced.
    Cell 200."""
    eps_grid = np.logspace(np.log10(eps_bounds[0]), np.log10(eps_bounds[1]), n_grid)
    mu2_grid = np.logspace(np.log10(mu2_bounds[0]), np.log10(mu2_bounds[1]), n_grid)
    frac_grid = np.logspace(np.log10(frac_bounds[0]), np.log10(frac_bounds[1]), n_grid)

    best = (-np.inf, eps_grid[0], mu2_grid[0], frac_grid[0])
    for e in eps_grid:
        for m2 in mu2_grid:
            for f in frac_grid:
                R = float(sns_key_rate_eta_joint(e, m2, f*m2, eta_value, e_a, p_dark=p_dark))
                if R > best[0]:
                    best = (R, e, m2, f)

    def neg(x):
        e, m2, f = x
        if not (eps_bounds[0] <= e <= eps_bounds[1]) or not (mu2_bounds[0] <= m2 <= mu2_bounds[1]) \
           or not (frac_bounds[0] <= f <= frac_bounds[1]):
            return 1e3
        return -float(sns_key_rate_eta_joint(e, m2, f*m2, eta_value, e_a, p_dark=p_dark))

    res = minimize(neg, x0=[best[1], best[2], best[3]], method="Nelder-Mead",
                    options={"xatol": 1e-12, "fatol": 1e-16, "maxiter": 3000, "maxfev": 3000})
    if res.success and -res.fun > best[0]:
        e, m2, f = res.x
        return float(-res.fun), float(e), float(m2), float(f*m2)
    return best[0], best[1], best[2], best[3]*best[2]


def optimize_joint_ensemble(eta_det, e_a, sigma_R2, p_dark, n_mc=20_000, seed=0, n_grid=10,
                             eps_bounds=(1e-5, 0.5), mu2_bounds=(1e-3, 1.0), frac_bounds=(0.005, 0.95)):
    """Ensemble (turbulence-faded) version of optimize_joint_signal_decoy.
    Reuses gamma_gamma_params/sample_gamma_gamma unchanged from Part II.
    Cell 202."""
    rng = np.random.default_rng(seed)
    alpha_g, beta_g = gamma_gamma_params(sigma_R2)
    h = sample_gamma_gamma(alpha_g, beta_g, n_mc, rng)
    eta_samples = eta_det * h

    eps_grid = np.logspace(np.log10(eps_bounds[0]), np.log10(eps_bounds[1]), n_grid)
    mu2_grid = np.logspace(np.log10(mu2_bounds[0]), np.log10(mu2_bounds[1]), n_grid)
    frac_grid = np.logspace(np.log10(frac_bounds[0]), np.log10(frac_bounds[1]), n_grid)
    best = (-np.inf, eps_grid[0], mu2_grid[0], frac_grid[0])
    for e in eps_grid:
        for m2 in mu2_grid:
            for f in frac_grid:
                R_mean = float(np.mean(sns_key_rate_eta_joint(e, m2, f*m2, eta_samples, e_a, p_dark=p_dark)))
                if R_mean > best[0]:
                    best = (R_mean, e, m2, f)

    def neg(x):
        e, m2, f = x
        if not (eps_bounds[0] <= e <= eps_bounds[1]) or not (mu2_bounds[0] <= m2 <= mu2_bounds[1]) \
           or not (frac_bounds[0] <= f <= frac_bounds[1]):
            return 1e3
        return -float(np.mean(sns_key_rate_eta_joint(e, m2, f*m2, eta_samples, e_a, p_dark=p_dark)))

    res = minimize(neg, x0=[best[1], best[2], best[3]], method="Nelder-Mead",
                    options={"xatol": 1e-10, "fatol": 1e-14, "maxiter": 1500})
    if res.success and -res.fun > best[0]:
        e, m2, f = res.x
        return float(-res.fun), float(e), float(m2), float(f*m2)
    return best[0], best[1], best[2], best[3]*best[2]
