"""Parts 5-7 — Optimisation.

* Part 5 (signal intensity mu'): instantaneous (perfect-CSI, per fading
  percentile) vs. ensemble (fixed operating point, Part 7's joint optimum).
* Part 6 (sending probability eps): a dedicated 1-D fine grid search holding
  mu' fixed, used for the eps-improvement-factor figure.
* Part 7 (joint (eps, mu')): the primary optimiser. Runs grid search
  (seed/validation), Differential Evolution (primary global search), and
  Basin Hopping (cross-check) against the Gauss-Laguerre ensemble-rate
  objective, and keeps the best of the three - exactly the triangulated
  approach requested in the task brief. The final reported rate is
  re-evaluated with a higher-order (80-node) quadrature for reporting
  precision.
"""

import numpy as np
from scipy.optimize import differential_evolution, basinhopping, minimize

from src.engine import sns_key_rate_eta
from src.fso.config import EPS_BOUNDS, MU_BOUNDS
from src.fso.channel_models import fso_deterministic_eta, gamma_gamma_params, sample_gamma_gamma
from src.fso.quadrature import gauss_laguerre_ensemble_rate


# ---------------------------------------------------------------------------
# Deterministic / instantaneous-channel optimisation (single fixed eta)
# ---------------------------------------------------------------------------
def optimize_instantaneous(eta_value, e_a, eps_bounds=EPS_BOUNDS, mu_bounds=MU_BOUNDS,
                            n_grid=20, seed=0):
    """Joint (eps, mu') optimisation of the key rate for ONE fixed, known
    channel transmittance eta_value (i.e. perfect-CSI / instantaneous-
    realization case). Coarse grid seed + Nelder-Mead polish (cheap: single
    scalar eta per call).
    """
    eps_grid = np.linspace(eps_bounds[0], eps_bounds[1], n_grid)
    mu_grid = np.linspace(mu_bounds[0], mu_bounds[1], n_grid)
    EPS, MU = np.meshgrid(eps_grid, mu_grid, indexing="ij")
    R = sns_key_rate_eta(EPS, MU, eta_value, e_a)
    idx = np.unravel_index(np.argmax(R), R.shape)
    best_eps, best_mu, best_R = EPS[idx], MU[idx], R[idx]

    def neg_rate(x):
        e, m = x
        if not (eps_bounds[0] <= e <= eps_bounds[1]) or not (mu_bounds[0] <= m <= mu_bounds[1]):
            return 1e3
        return -float(sns_key_rate_eta(e, m, eta_value, e_a))

    res = minimize(neg_rate, x0=[best_eps, best_mu], method="Nelder-Mead",
                    options={"xatol": 1e-10, "fatol": 1e-14, "maxiter": 600})
    if res.success and -res.fun > best_R:
        return float(-res.fun), float(res.x[0]), float(res.x[1])
    return float(best_R), float(best_eps), float(best_mu)


def instantaneous_maps(distances_km, e_a, sigma_R2_fn, alpha_fso=None, eta_d=0.8,
                        percentiles=(5, 25, 50, 75, 95), n_mc_percentile=200_000, seed=0):
    """Part 5/6 deliverable: 'instantaneous alpha/epsilon map'.
    For every distance and every percentile of the (distance-dependent)
    Gamma-Gamma fading distribution, find the instantaneous eta realization
    at that percentile and jointly-optimise (eps, mu') for it with perfect
    CSI.

    sigma_R2_fn: callable d_km -> sigma_R^2 (so the caller controls which
    turbulence case, e.g. weak/moderate/strong, is used to build the fading
    distribution).

    Returns a DataFrame-ready list of dict rows.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for d in distances_km:
        eta_det = fso_deterministic_eta(d, eta_d=eta_d) if alpha_fso is None \
            else fso_deterministic_eta(d, alpha=alpha_fso, eta_d=eta_d)
        sigma2 = sigma_R2_fn(d)
        alpha_g, beta_g = gamma_gamma_params(sigma2)
        h_samples = sample_gamma_gamma(alpha_g, beta_g, n_mc_percentile, rng)
        for p in percentiles:
            h_p = float(np.percentile(h_samples, p))
            eta_inst = eta_det * h_p
            R_opt, eps_opt, mu_opt = optimize_instantaneous(eta_inst, e_a)
            rows.append(dict(d_km=d, percentile=p, h_p=h_p, eta_inst=eta_inst,
                              R_inst=R_opt, eps_inst=eps_opt, mu_inst=mu_opt))
    return rows


# ---------------------------------------------------------------------------
# Ensemble-averaged optimisation (Part 7: joint eps,mu' over Gamma-Gamma channel)
# ---------------------------------------------------------------------------
def _neg_ensemble_rate(x, d_km, e_a, sigma_R2, n_points, alpha_fso, eta_d):
    eps, mu = x
    if not (EPS_BOUNDS[0] <= eps <= EPS_BOUNDS[1]) or not (MU_BOUNDS[0] <= mu <= MU_BOUNDS[1]):
        return 1e3
    return -gauss_laguerre_ensemble_rate(eps, mu, d_km, e_a, sigma_R2, n_points, alpha_fso, eta_d)


def grid_search_joint_ensemble(d_km, e_a, sigma_R2, n_grid=20, n_points=20,
                                alpha_fso=None, eta_d=0.8):
    """Part 7 validation method: coarse grid search over (eps, mu')."""
    eps_grid = np.linspace(*EPS_BOUNDS, n_grid)
    mu_grid = np.linspace(*MU_BOUNDS, n_grid)
    best = (-np.inf, eps_grid[0], mu_grid[0])
    for e in eps_grid:
        for m in mu_grid:
            R = gauss_laguerre_ensemble_rate(e, m, d_km, e_a, sigma_R2, n_points, alpha_fso, eta_d)
            if R > best[0]:
                best = (R, e, m)
    return best  # (R, eps, mu)


def differential_evolution_joint_ensemble(d_km, e_a, sigma_R2, n_points=20,
                                           alpha_fso=None, eta_d=0.8, seed=0, record_history=False):
    """Part 7 primary global optimiser: Differential Evolution."""
    history = []

    def cb(xk, convergence):
        if record_history:
            R = gauss_laguerre_ensemble_rate(xk[0], xk[1], d_km, e_a, sigma_R2, n_points, alpha_fso, eta_d)
            history.append((xk[0], xk[1], R))

    res = differential_evolution(
        _neg_ensemble_rate, bounds=[EPS_BOUNDS, MU_BOUNDS],
        args=(d_km, e_a, sigma_R2, n_points, alpha_fso, eta_d),
        seed=seed, maxiter=150, popsize=20, tol=1e-10, polish=True,
        callback=cb if record_history else None, updating="deferred",
    )
    out = (-res.fun, res.x[0], res.x[1])
    if record_history:
        return out, history
    return out


def basin_hopping_joint_ensemble(d_km, e_a, sigma_R2, x0, n_points=20,
                                  alpha_fso=None, eta_d=0.8, seed=0):
    """Part 7 cross-check global optimiser: Basin Hopping (Nelder-Mead local
    steps)."""
    minimizer_kwargs = dict(method="Nelder-Mead",
                             args=(d_km, e_a, sigma_R2, n_points, alpha_fso, eta_d),
                             options={"xatol": 1e-9, "fatol": 1e-13})
    res = basinhopping(_neg_ensemble_rate, x0=x0, minimizer_kwargs=minimizer_kwargs,
                        niter=40, seed=seed, stepsize=0.15)
    eps_opt = float(np.clip(res.x[0], *EPS_BOUNDS))
    mu_opt = float(np.clip(res.x[1], *MU_BOUNDS))
    R_opt = gauss_laguerre_ensemble_rate(eps_opt, mu_opt, d_km, e_a, sigma_R2, n_points, alpha_fso, eta_d)
    return (R_opt, eps_opt, mu_opt)


def joint_optimize_ensemble(d_km, e_a, sigma_R2, n_points=20, alpha_fso=None, eta_d=0.8,
                             seed=0, record_history=False):
    """Part 7 top-level driver: run grid search (seed/validation),
    Differential Evolution (primary global search), and Basin Hopping
    (cross-check), then return the best of the three together with a small
    comparison table - exactly the triangulated-optimum approach requested
    in the task brief. Final reported rate is re-evaluated with a
    higher-order (n_points=40) quadrature for reporting precision.
    """
    grid_best = grid_search_joint_ensemble(d_km, e_a, sigma_R2, n_grid=20, n_points=n_points,
                                            alpha_fso=alpha_fso, eta_d=eta_d)
    de_result = differential_evolution_joint_ensemble(d_km, e_a, sigma_R2, n_points=n_points,
                                                       alpha_fso=alpha_fso, eta_d=eta_d, seed=seed,
                                                       record_history=record_history)
    de_best, history = (de_result if record_history else (de_result, None))
    bh_best = basin_hopping_joint_ensemble(d_km, e_a, sigma_R2, x0=[de_best[1], de_best[2]],
                                            n_points=n_points, alpha_fso=alpha_fso, eta_d=eta_d, seed=seed)

    candidates = {"grid": grid_best, "differential_evolution": de_best, "basin_hopping": bh_best}
    best_method = max(candidates, key=lambda k: candidates[k][0])
    R_best, eps_best, mu_best = candidates[best_method]

    # Precision re-evaluation at n_points=40
    R_final = gauss_laguerre_ensemble_rate(eps_best, mu_best, d_km, e_a, sigma_R2, 40, alpha_fso, eta_d)

    out = dict(d_km=d_km, R=R_final, eps=eps_best, mu=mu_best, best_method=best_method,
               R_grid=grid_best[0], eps_grid=grid_best[1], mu_grid=grid_best[2],
               R_de=de_best[0], eps_de=de_best[1], mu_de=de_best[2],
               R_bh=bh_best[0], eps_bh=bh_best[1], mu_bh=bh_best[2])
    if record_history:
        out["history"] = history
    return out


# ---------------------------------------------------------------------------
# Part 6: epsilon-only optimisation (mu' held at the ensemble joint optimum)
# ---------------------------------------------------------------------------
def optimize_epsilon_ensemble(d_km, mu_fixed, e_a, sigma_R2, n_points=40,
                               alpha_fso=None, eta_d=0.8, n_grid=200):
    """Fine 1-D grid search for eps* holding mu' fixed (Part 6)."""
    eps_grid = np.linspace(*EPS_BOUNDS, n_grid)
    R_grid = np.array([gauss_laguerre_ensemble_rate(e, mu_fixed, d_km, e_a, sigma_R2, n_points, alpha_fso, eta_d)
                        for e in eps_grid])
    i = int(np.argmax(R_grid))
    return float(R_grid[i]), float(eps_grid[i])
