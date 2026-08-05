"""Part 4 — Gauss-Laguerre Quadrature Validation.

Since h = x*y/(alpha_g*beta_g) with x~Gamma(alpha_g,1), y~Gamma(beta_g,1), the
ensemble average <R> is exactly a double generalized Gauss-Laguerre
quadrature. Deep in the saturation regime alpha_g can reach the hundreds,
where `roots_genlaguerre`'s raw weights overflow float64 (they scale like
Gamma(alpha_g)); this is handled with (i) per-factor weight renormalisation
instead of dividing by exp(gammaln(.)), and (ii) a Gauss-Hermite fallback
using the Gaussian (CLT) limit Gamma(k,1) ~ N(k,k) for shape parameters above
~150.
"""

import numpy as np
from scipy.special import roots_genlaguerre

from src.engine import sns_key_rate_eta, full_metrics_eta_vec
from src.fso.config import N_MC_DEFAULT
from src.fso.channel_models import fso_deterministic_eta, gamma_gamma_params
from src.fso.monte_carlo import monte_carlo_ensemble_rate

_GL_CACHE = {}

# scipy.special.roots_genlaguerre becomes numerically unstable (weights overflow
# float64, since they scale like Gamma(shape)) once shape exceeds ~150-180. Deep
# in the turbulence saturation regime, alpha_g can reach the hundreds (see Part 2),
# so for shape parameters above GL_SHAPE_OVERFLOW_LIMIT we instead use the
# Gaussian (CLT) limit of Gamma(k,1) = k + sqrt(k)*Z, Z~N(0,1), evaluated at
# Gauss-Hermite nodes. This is an extremely accurate approximation in this
# regime (relative skewness ~ 2/sqrt(k) < 0.15 for k>180) and is the standard
# fix for this well-known instability of generalized Gauss-Laguerre quadrature.
GL_SHAPE_OVERFLOW_LIMIT = 150.0


def _nodes_weights_for_gamma_factor(shape, n_points):
    """1-D nodes/weights (summing to 1) approximating E_X[g(X)] for
    X~Gamma(shape,1)."""
    if shape < GL_SHAPE_OVERFLOW_LIMIT:
        x, w = roots_genlaguerre(n_points, shape - 1)
        w = w / w.sum()
        return x, w
    else:
        # Gaussian (CLT) limit: X ~ N(shape, shape)
        z, wz = np.polynomial.hermite_e.hermegauss(n_points)  # weight exp(-z^2/2)
        wz = wz / wz.sum()
        x = shape + np.sqrt(shape) * z
        x = np.clip(x, 1e-12, None)  # guard against (astronomically unlikely) negative tail
        return x, wz


def _gl_nodes_weights(sigma_R2, n_points):
    key = (round(sigma_R2, 10), n_points)
    if key not in _GL_CACHE:
        alpha_g, beta_g = gamma_gamma_params(sigma_R2)
        x, nwx = _nodes_weights_for_gamma_factor(alpha_g, n_points)
        y, nwy = _nodes_weights_for_gamma_factor(beta_g, n_points)
        X, Y = np.meshgrid(x, y, indexing="ij")
        NWX, NWY = np.meshgrid(nwx, nwy, indexing="ij")
        H = (X * Y / (alpha_g * beta_g)).ravel()
        W = (NWX * NWY).ravel()  # already normalised, sums to 1
        _GL_CACHE[key] = (H, W)
    return _GL_CACHE[key]


def gauss_laguerre_ensemble_rate(eps, mu_prime, d_km, e_a, sigma_R2, n_points=40,
                                  alpha_fso=None, eta_d=0.8):
    """Double generalized Gauss-Laguerre quadrature estimate of
        <R> = int_0^inf R(eta_det*h) f_GG(h) dh
    using the exact representation h = x*y/(alpha_g*beta_g), x~Gamma(alpha_g,1),
    y~Gamma(beta_g,1):

        <R> = 1/(Gamma(alpha_g)*Gamma(beta_g)) * sum_i sum_j w_i*w_j*R(eta_det*x_i*y_j/(alpha_g*beta_g))

    where {x_i,w_i} and {y_j,w_j} are the n_points-node generalized
    Gauss-Laguerre nodes/weights for weight functions x^(alpha_g-1)e^-x and
    y^(beta_g-1)e^-y respectively (scipy.special.roots_genlaguerre), each
    independently renormalised to sum to 1 rather than divided by
    exp(gammaln(.)) (see GL_SHAPE_OVERFLOW_LIMIT note above for why, and the
    Gauss-Hermite fallback used above that limit).
    """
    H, W = _gl_nodes_weights(sigma_R2, n_points)

    eta_det = fso_deterministic_eta(d_km, eta_d=eta_d) if alpha_fso is None \
        else fso_deterministic_eta(d_km, alpha=alpha_fso, eta_d=eta_d)
    R_nodes = sns_key_rate_eta(eps, mu_prime, eta_det * H, e_a)
    return float(np.sum(W * R_nodes))


def gauss_laguerre_ensemble_metrics(eps, mu_prime, d_km, e_a, sigma_R2, n_points=40,
                                     alpha_fso=None, eta_d=0.8):
    """Ensemble-averaged version of EVERY performance metric (Part 8), not
    just R, via the same double Gauss-Laguerre quadrature used for the key
    rate."""
    H, W = _gl_nodes_weights(sigma_R2, n_points)
    eta_det = fso_deterministic_eta(d_km, eta_d=eta_d) if alpha_fso is None \
        else fso_deterministic_eta(d_km, alpha=alpha_fso, eta_d=eta_d)
    metrics = full_metrics_eta_vec(eps, mu_prime, eta_det * H, e_a)
    return {k: float(np.sum(W * v)) for k, v in metrics.items()}


def quadrature_validation(eps, mu_prime, d_km, e_a, sigma_R2, n_mc=N_MC_DEFAULT,
                           n_points_list=(20, 40, 80), alpha_fso=None, eta_d=0.8, seed=0):
    """Part 4 deliverable: compare MC estimate against 20/40/80-point GL
    quadrature."""
    rng = np.random.default_rng(seed)
    mc = monte_carlo_ensemble_rate(eps, mu_prime, d_km, e_a, sigma_R2, n_mc=n_mc,
                                    alpha_fso=alpha_fso, eta_d=eta_d, rng=rng)
    rows = []
    for n in n_points_list:
        R_gl = gauss_laguerre_ensemble_rate(eps, mu_prime, d_km, e_a, sigma_R2, n_points=n,
                                             alpha_fso=alpha_fso, eta_d=eta_d)
        abs_err = abs(mc["mean"] - R_gl)
        rel_err = abs_err / abs(mc["mean"]) if mc["mean"] != 0 else np.nan
        rows.append(dict(n_points=n, R_GL=R_gl, R_MC=mc["mean"], abs_err=abs_err, rel_err=rel_err))
    return mc, rows
