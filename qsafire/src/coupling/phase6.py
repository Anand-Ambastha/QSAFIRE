"""Coupling to the SNS-TF-QKD engine (notebook Part VIII / Phase 6).

Folds Phase 4/5's turbulence fading statistics (Gamma-Gamma) properly into
the Phase 5 real-loss deterministic channel via Monte Carlo ensemble
optimisation, with the real (suppressed) p_dark threaded through explicitly
-- mirrors Part II's `joint_optimize_ensemble` logic but exposes p_dark as a
pass-through parameter, which the original does not.

Origin: notebook cell 167 (optimize_ensemble_real).
"""

import numpy as np
from scipy.optimize import minimize

from src.engine import sns_key_rate_eta
from src.fso.channel_models import gamma_gamma_params, sample_gamma_gamma


def optimize_ensemble_real(eta_det, e_a, sigma_R2, p_dark, n_mc=60_000, seed=0,
                            eps_bounds=(0.01, 0.99), mu_bounds=(0.01, 1.0), n_grid=14):
    """Ensemble (turbulence-faded) joint (eps, mu') optimisation, real p_dark
    threaded through. Reuses gamma_gamma_params / sample_gamma_gamma
    unchanged from Part II. Cell 167."""
    rng = np.random.default_rng(seed)
    alpha_g, beta_g = gamma_gamma_params(sigma_R2)
    h = sample_gamma_gamma(alpha_g, beta_g, n_mc, rng)
    eta_samples = eta_det * h

    eps_grid = np.linspace(*eps_bounds, n_grid)
    mu_grid = np.linspace(*mu_bounds, n_grid)
    best = (-np.inf, eps_grid[0], mu_grid[0])
    for e in eps_grid:
        for m in mu_grid:
            R_mean = float(np.mean(sns_key_rate_eta(e, m, eta_samples, e_a, p_dark=p_dark)))
            if R_mean > best[0]:
                best = (R_mean, e, m)

    def neg_rate(x):
        e, m = x
        if not (eps_bounds[0] <= e <= eps_bounds[1]) or not (mu_bounds[0] <= m <= mu_bounds[1]):
            return 1e3
        return -float(np.mean(sns_key_rate_eta(e, m, eta_samples, e_a, p_dark=p_dark)))

    res = minimize(neg_rate, x0=[best[1], best[2]], method="Nelder-Mead",
                    options={"xatol": 1e-9, "fatol": 1e-13, "maxiter": 500})
    if res.success and -res.fun > best[0]:
        return float(-res.fun), float(res.x[0]), float(res.x[1])
    return best
