"""Fully asymmetric SNS-TF-QKD rate formulas, derived from first principles
(notebook Part XIII / Phase 9 final).

Two wrong approaches were tried and rejected before this: (1) per-station
rates -- wrong, because SNS-TF-QKD is a joint three-party protocol with one
shared Charlie measurement, not two independent point-to-point links; (2) a
straight product eta_A*eta_B as a single effective eta -- wrong, because
twin-field protocols get their sqrt-of-total-loss scaling advantage
precisely from each arm contributing *linearly* (not multiplicatively) to
Charlie's click probability; multiplying first would silently erase that
advantage. The formulas below keep each arm's transmittance separate all the
way through, splitting the Z-window into its four true cases (Alice-only,
Bob-only, both, neither).

Origin: notebook cell 240 (gain_X_asym, error_X_asym, S_Z_gain_asym,
E_Z_qber_asym, sns_key_rate_asym_full + symmetric-limit validation), cell 242
(CITY_BG_MULTIPLIER, optimize_joint_asym_ensemble).
"""

import numpy as np
from scipy.optimize import minimize

from src.engine import (
    P_DARK, F_EC, s0_vacuum_yield, s1_decoy_lower_bound, e1ph_upper_bound, binary_entropy,
    gain_X, S_Z_gain, E_Z_qber,
)
from src.fso.channel_models import gamma_gamma_params, sample_gamma_gamma

CITY_BG_MULTIPLIER = {
    "Delhi":     dict(mult=2.0, conf="Medium (direction real, Bedi et al. 2021; magnitude Indicative)"),
    "Mumbai":    dict(mult=1.5, conf="Medium (direction real, Bedi et al. 2021; magnitude Indicative)"),
    "Nainital":  dict(mult=0.3, conf="Indicative (small town / observatory site, no measurement)"),
    "Ahmedabad": dict(mult=1.5, conf="Indicative (generic large-city default, no measurement)"),
    "Bangalore": dict(mult=1.5, conf="Indicative (generic large-city default, no measurement)"),
    "Hyderabad": dict(mult=1.5, conf="Indicative (generic large-city default, no measurement)"),
}


def gain_X_asym(mu, eta_A, eta_B, p_dark=P_DARK):
    """Cell 240."""
    lam = mu*eta_A + mu*eta_B
    return 1.0 - (1.0 - 2.0*p_dark)*np.exp(-lam)


def error_X_asym(mu, eta_A, eta_B, e_a, p_dark=P_DARK):
    """Cell 240."""
    lam = mu*eta_A + mu*eta_B
    correct_click_prob = 1.0 - np.exp(-lam)
    dark_click_prob = 2.0*p_dark*np.exp(-lam)
    S = gain_X_asym(mu, eta_A, eta_B, p_dark)
    numerator = e_a*correct_click_prob + 0.5*dark_click_prob
    return np.where(S > 0, numerator/np.where(S > 0, S, 1.0), 0.5)


def S_Z_gain_asym(eps, mu2, eta_A, eta_B, p_dark=P_DARK):
    """Cell 240."""
    alice_only_click = 1.0 - np.exp(-mu2*eta_A)
    bob_only_click    = 1.0 - np.exp(-mu2*eta_B)
    both_send_click   = 1.0 - np.exp(-(mu2*eta_A + mu2*eta_B))
    neither_send_click = 2.0*p_dark
    return (eps*(1-eps)*alice_only_click + (1-eps)*eps*bob_only_click
            + eps**2*both_send_click + (1-eps)**2*neither_send_click)


def E_Z_qber_asym(eps, mu2, eta_A, eta_B, p_dark=P_DARK):
    """Cell 240."""
    both_send_click = 1.0 - np.exp(-(mu2*eta_A + mu2*eta_B))
    neither_send_click = 2.0*p_dark
    Sz = S_Z_gain_asym(eps, mu2, eta_A, eta_B, p_dark)
    numerator = eps**2*both_send_click + (1-eps)**2*neither_send_click
    return np.where(Sz > 0, numerator/np.where(Sz > 0, Sz, 1.0), 0.5)


def sns_key_rate_asym_full(eps, mu2, mu1, eta_A, eta_B, e_a, f_ec=F_EC, p_dark=P_DARK):
    """Cell 240."""
    eta_A = np.asarray(eta_A, dtype=float); eta_B = np.asarray(eta_B, dtype=float)
    s0 = s0_vacuum_yield(p_dark)
    S1 = gain_X_asym(mu1, eta_A, eta_B, p_dark)
    S2 = gain_X_asym(mu2, eta_A, eta_B, p_dark)
    E1 = error_X_asym(mu1, eta_A, eta_B, e_a, p_dark)
    s1 = np.maximum(s1_decoy_lower_bound(mu1, mu2, S1, S2, s0), 1e-15)
    e1ph = np.clip(e1ph_upper_bound(mu1, S1, E1, s0, s1), 0.0, 0.5)
    Sz = S_Z_gain_asym(eps, mu2, eta_A, eta_B, p_dark)
    Ez = np.clip(E_Z_qber_asym(eps, mu2, eta_A, eta_B, p_dark), 1e-12, 0.5)
    R = (2*eps*(1-eps)*mu2*np.exp(-mu2)*s1*(1-binary_entropy(e1ph)) - Sz*f_ec*binary_entropy(Ez))
    return dict(S1=S1, S2=S2, E1=E1, s0=s0, s1=s1, e1ph=e1ph, Sz=Sz, Ez=Ez, R=R)


def validate_symmetric_limit(eta_t=0.01, p_dark=P_DARK):
    """Cell 240 (script part): setting eta_A = eta_B must reproduce the
    original symmetric Part I formulas to machine precision."""
    checks = []
    for mu2_t, mu1_t, eps_t in [(0.3, 0.06, 0.05), (0.41, 0.002, 0.018), (0.1, 0.02, 0.1)]:
        g_old, g_new = gain_X(mu1_t, eta_t, p_dark), gain_X_asym(mu1_t, eta_t, eta_t, p_dark)
        sz_old, sz_new = S_Z_gain(eps_t, mu2_t, eta_t, p_dark), S_Z_gain_asym(eps_t, mu2_t, eta_t, eta_t, p_dark)
        ez_old, ez_new = E_Z_qber(eps_t, mu2_t, eta_t, p_dark), E_Z_qber_asym(eps_t, mu2_t, eta_t, eta_t, p_dark)
        checks.append(max(abs(g_old-g_new), abs(sz_old-sz_new), abs(ez_old-ez_new)))
    max_dev = max(checks)
    print(f"Max deviation from original symmetric formulas across 3 test points: {max_dev:.2e}")
    assert max_dev < 1e-12, "Asymmetric formulas do not reduce to the symmetric case exactly"
    print("PASS: exact reduction confirmed.")
    return max_dev


def optimize_joint_asym_ensemble(eta_A_det, eta_B_det, sigma_A, sigma_B, e_a, p_dark,
                                  n_mc=4000, n_grid=6, seed=0,
                                  eps_bounds=(1e-5, 0.5), mu2_bounds=(1e-3, 1.0), frac_bounds=(0.005, 0.95)):
    """Cell 242: independent per-arm Gamma-Gamma fading (Delhi and Mumbai,
    for example, have different Rytov variances and should not share one
    averaged fading process)."""
    rng_A = np.random.default_rng(seed)
    rng_B = np.random.default_rng(seed + 1)  # independent stream -- arms fade independently
    agA, bgA = gamma_gamma_params(sigma_A)
    agB, bgB = gamma_gamma_params(sigma_B)
    hA = sample_gamma_gamma(agA, bgA, n_mc, rng_A)
    hB = sample_gamma_gamma(agB, bgB, n_mc, rng_B)
    etaA_s, etaB_s = eta_A_det*hA, eta_B_det*hB

    eps_grid = np.logspace(np.log10(eps_bounds[0]), np.log10(eps_bounds[1]), n_grid)
    mu2_grid = np.logspace(np.log10(mu2_bounds[0]), np.log10(mu2_bounds[1]), n_grid)
    frac_grid = np.logspace(np.log10(frac_bounds[0]), np.log10(frac_bounds[1]), n_grid)
    best = (-np.inf, eps_grid[0], mu2_grid[0], frac_grid[0])
    for e in eps_grid:
        for m2 in mu2_grid:
            for f in frac_grid:
                Rm = float(np.mean(sns_key_rate_asym_full(e, m2, f*m2, etaA_s, etaB_s, e_a, p_dark=p_dark)['R']))
                if Rm > best[0]:
                    best = (Rm, e, m2, f)

    def neg(x):
        e, m2, f = x
        if not (eps_bounds[0] <= e <= eps_bounds[1]) or not (mu2_bounds[0] <= m2 <= mu2_bounds[1]) \
           or not (frac_bounds[0] <= f <= frac_bounds[1]):
            return 1e3
        return -float(np.mean(sns_key_rate_asym_full(e, m2, f*m2, etaA_s, etaB_s, e_a, p_dark=p_dark)['R']))

    res = minimize(neg, x0=[best[1], best[2], best[3]], method="Nelder-Mead",
                    options={"xatol": 1e-9, "fatol": 1e-13, "maxiter": 800})
    if res.success and -res.fun > best[0]:
        e, m2, f = res.x
        return float(-res.fun), float(e), float(m2), float(f*m2)
    return best[0], best[1], best[2], best[3]*best[2]
