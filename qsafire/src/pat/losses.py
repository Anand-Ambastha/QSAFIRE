"""PAT, optics and detector real-loss module (notebook Part VII / Phase 5).

Replaces the idealized ETA_D=0.8 / P_DARK=1e-11 placeholders used throughout
Parts I-VI with real, individually-cited hardware numbers: Friis-equation
diffraction/collection efficiency, Gaussian-beam pointing loss, a real
characterised Si APD, and day/night + city-specific background-light
scenarios.

Origin: notebook cells 145 (hardware baseline + eta_diffraction), 147
(eta_pointing, POINTING_CASES), 149 (detector realism constants), 151
(R_NOISE_BASELINE_HZ, BACKGROUND_SCENARIOS), 153
(optimize_instantaneous_real).
"""

import numpy as np
from scipy.optimize import minimize

from src.engine import sns_key_rate_eta

# ============================================================
# 7.2 Hardware baseline (Micius-type transmitter + 1 m ground receiver)
# ============================================================
LAMBDA_P5_M = 810e-9          # match Phases 2-4
THETA_DIV_RAD = 10e-6         # satellite Tx half-divergence (Liao et al. 2017, Phase 1 doc Section 5.1)
D_RX_M = 1.0                  # ground receiver primary (Micius-style, Phase 1 doc Section 5.1)
OBSCURATION_RATIO = 0.3/0.8   # DLR OGS-KN secondary/primary ratio [G2012], applied to our own aperture
D_RX_INT_M = D_RX_M * OBSCURATION_RATIO
R_TX_HZ = 100e6               # assumed system pulse rate (Micius-representative; order-of-magnitude, Indicative)

G_TX = 8.0 / THETA_DIV_RAD**2
G_RX = (np.pi**2 / LAMBDA_P5_M**2) * (D_RX_M**2 - D_RX_INT_M**2)

ETA_RX_INTERNAL = 0.95**7     # 7 mirrors at 95% reflectivity each [S2025]

POINTING_CASES = {
    "Measured (Micius, [W2021])":    dict(sigma_j=0.47e-6, conf="High (real in-orbit measurement)"),
    "Conservative default ([S2025])": dict(sigma_j=2.5e-6,  conf="Medium (literature conservative design value)"),
}

# ============================================================
# 7.4 Detector realism
# ============================================================
ETA_DET_REAL = 0.50           # photon detection efficiency, Si APD [A2021]
DARK_COUNT_RATE_HZ = 1.0      # conservative upper bound on "< 1 Hz" [A2021]
TIMING_JITTER_NS = 1.0        # < 1 ns [A2021]
AFTERPULSE_PROB = 0.0035      # < 0.35% [A2021]
P_DARK_REAL = DARK_COUNT_RATE_HZ / R_TX_HZ

# ============================================================
# 7.5 Background light
# ============================================================
R_NOISE_BASELINE_HZ = 3000.0   # representative "typical working" background rate, [S2025] Table I

BACKGROUND_SCENARIOS = {
    "Delhi (Alice) -- night": dict(mult=2.0,    conf="Indicative (direction real, [B2021]; magnitude estimated)"),
    "Mumbai (Bob) -- night":  dict(mult=1.5,    conf="Indicative (direction real, [B2021]; magnitude estimated)"),
    "Either city -- day":     dict(mult=1000.0, conf="Indicative (order-of-magnitude qualitative only)"),
}

E_A_MISALIGNMENT = 0.15   # retained from the existing notebook's default misalignment sweep


def eta_diffraction(L_m, G_tx=G_TX, G_rx=G_RX, lam=LAMBDA_P5_M):
    """Friis-equation collection efficiency [S2025], Eq. (12)-(13). Cell 145."""
    L_freespace = (lam / (4 * np.pi * L_m))**2
    return G_tx * G_rx * L_freespace


def eta_pointing(sigma_jitter_rad, theta_div=THETA_DIV_RAD, theta_bias=0.0):
    """Gaussian-beam pointing-loss formula ([S2025] Eq. (14), citing [CF2023]). Cell 147."""
    denom = theta_div**2 + 4*sigma_jitter_rad**2
    return (theta_div**2/denom) * np.exp(-2*theta_bias**2/denom)


def optimize_instantaneous_real(eta_value, e_a, p_dark, eps_bounds=(0.01, 0.99), mu_bounds=(0.01, 1.0), n_grid=20):
    """Local optimizer that threads a real p_dark override through
    sns_key_rate_eta without relying on the module-level default (which is
    bound at function-definition time). Cell 153."""
    eps_grid = np.linspace(*eps_bounds, n_grid)
    mu_grid = np.linspace(*mu_bounds, n_grid)
    EPS, MU = np.meshgrid(eps_grid, mu_grid, indexing="ij")
    R = sns_key_rate_eta(EPS, MU, eta_value, e_a, p_dark=p_dark)
    idx = np.unravel_index(np.argmax(R), R.shape)
    best_eps, best_mu, best_R = EPS[idx], MU[idx], R[idx]

    def neg_rate(x):
        e, m = x
        if not (eps_bounds[0] <= e <= eps_bounds[1]) or not (mu_bounds[0] <= m <= mu_bounds[1]):
            return 1e3
        return -float(sns_key_rate_eta(e, m, eta_value, e_a, p_dark=p_dark))

    res = minimize(neg_rate, x0=[best_eps, best_mu], method="Nelder-Mead",
                    options={"xatol": 1e-10, "fatol": 1e-14, "maxiter": 600})
    if res.success and -res.fun > best_R:
        return float(-res.fun), float(res.x[0]), float(res.x[1])
    return float(best_R), float(best_eps), float(best_mu)
