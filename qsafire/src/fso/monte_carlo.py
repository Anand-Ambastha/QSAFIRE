"""Part 3 — Monte Carlo Turbulence Propagation.

Draws N_MC Gamma-Gamma fading samples, forms the effective per-realization
transmittance, and evaluates the SNS-TF-QKD key rate for every realization
(fully vectorised, no Python loop over samples).
"""

import numpy as np

from src.engine import sns_key_rate_eta
from src.fso.config import N_MC_DEFAULT
from src.fso.channel_models import fso_deterministic_eta, gamma_gamma_params, sample_gamma_gamma


def monte_carlo_ensemble_rate(eps, mu_prime, d_km, e_a, sigma_R2, n_mc=N_MC_DEFAULT,
                               alpha_fso=None, eta_d=0.8, rng=None, return_samples=False):
    """Draw n_mc Gamma-Gamma fading samples h_i, form the effective
    per-realization transmittance eta_i = eta_det(d) * h_i, evaluate
    R_i = R(eps, mu', eta_i, e_a) for EVERY realization (vectorised - no
    Python loop over samples), and return Monte Carlo summary statistics of
    R.
    """
    if rng is None:
        rng = np.random.default_rng()
    eta_det = fso_deterministic_eta(d_km, eta_d=eta_d) if alpha_fso is None \
        else fso_deterministic_eta(d_km, alpha=alpha_fso, eta_d=eta_d)
    alpha_g, beta_g = gamma_gamma_params(sigma_R2)
    h = sample_gamma_gamma(alpha_g, beta_g, n_mc, rng)
    eta_i = eta_det * h
    R_i = sns_key_rate_eta(eps, mu_prime, eta_i, e_a)

    mean_R = float(np.mean(R_i))
    std_R = float(np.std(R_i, ddof=1))
    sem = std_R / np.sqrt(n_mc)
    ci95 = (mean_R - 1.96 * sem, mean_R + 1.96 * sem)

    out = dict(mean=mean_R, var=float(np.var(R_i, ddof=1)), std=std_R,
               sem=sem, ci95_lo=ci95[0], ci95_hi=ci95[1], n_mc=n_mc,
               eta_det=float(eta_det))
    if return_samples:
        out["R_samples"] = R_i
        out["h_samples"] = h
    return out


def monte_carlo_convergence(eps, mu_prime, d_km, e_a, sigma_R2, n_list, alpha_fso=None,
                             eta_d=0.8, seed=0):
    """Part 9 diagnostic: running-mean convergence of the MC estimator vs
    N_MC."""
    rng = np.random.default_rng(seed)
    eta_det = fso_deterministic_eta(d_km, eta_d=eta_d) if alpha_fso is None \
        else fso_deterministic_eta(d_km, alpha=alpha_fso, eta_d=eta_d)
    alpha_g, beta_g = gamma_gamma_params(sigma_R2)
    n_max = max(n_list)
    h_full = sample_gamma_gamma(alpha_g, beta_g, n_max, rng)
    R_full = sns_key_rate_eta(eps, mu_prime, eta_det * h_full, e_a)
    running_mean = np.cumsum(R_full) / np.arange(1, n_max + 1)
    return np.array(n_list), running_mean[np.array(n_list) - 1]
