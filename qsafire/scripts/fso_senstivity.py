"""One-at-a-time (OAT) sensitivity analysis for the Part I/II FSO link
(src/fso, src/engine).

Two separate sweeps, because the two families of parameters enter the model
through different code paths:

(A) "Deterministic-channel" parameters (d_km, alpha_fso, eta_d, e_a,
    decoy_ratio, f_ec, p_dark): these determine eta_det and the engine's rate
    formula directly. Swept against a single fixed eta_det, jointly
    reoptimizing (eps, mu') at every perturbation via a local optimizer that
    (unlike optimize_instantaneous in src/fso/optimization.py) exposes
    decoy_ratio/f_ec/p_dark as free arguments, since those are only default
    kwargs on sns_key_rate_eta and are not threaded through the existing
    wrapper.

(B) The turbulence parameter Cn2: only enters through the Gamma-Gamma
    ensemble path (sigma_R2 -> alpha_g/beta_g -> h -> eta_i), so it is swept
    using the existing, already-verified joint_optimize_ensemble ensemble
    optimizer at the SAME baseline distance/alpha/eta_d/e_a as (A). Values
    of Cn2 or d_km that would push sigma_R2 past
    channel_models.SIGMA_R2_VALIDITY_MAX are skipped and flagged, per the
    validity-range fix in channel_models.py.

For every parameter, perturbations of -20%, -10%, +10%, +20% around the
baseline are evaluated one at a time (all other parameters held at
baseline), and the elasticity

    E = (%change in R) / (%change in parameter)

is reported as the headline sensitivity measure so that parameters with very
different natural units/scales (e.g. p_dark ~1e-11 vs eta_d ~0.8) are
directly comparable. Rows are sorted by |E| at +-10%, descending.

Run with:  python -m scripts.sensitivity_fso
Outputs:   outputs/csv/sensitivity_fso.csv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from src.engine import sns_key_rate_eta, DECOY_RATIO, F_EC, P_DARK
from src.fso.config import ALPHA_FSO_DB_KM, CN2_MODERATE, E_A
from src.fso.channel_models import (
    fso_deterministic_eta, rytov_variance, gamma_gamma_params,
    sigma_R2_within_validity, SIGMA_R2_VALIDITY_MAX,
)
from src.fso.optimization import joint_optimize_ensemble

OUT_CSV = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "outputs", "csv", "sensitivity_fso.csv")

# ---------------------------------------------------------------------------
# Baseline operating point
# ---------------------------------------------------------------------------
BASE = dict(
    d_km=5.0,               # within the GG model's validity range at CN2_MODERATE
    alpha_fso=ALPHA_FSO_DB_KM,
    eta_d=0.8,
    e_a=E_A,
    decoy_ratio=DECOY_RATIO,
    f_ec=F_EC,
    p_dark=P_DARK,
    Cn2=CN2_MODERATE,
)
PERTURBATIONS = (-0.20, -0.10, 0.10, 0.20)


def optimize_R_full(eta_value, e_a, decoy_ratio, f_ec, p_dark, n_grid=25):
    """Joint (eps, mu') optimizer for a single fixed eta that exposes
    decoy_ratio/f_ec/p_dark as free arguments (grid seed + Nelder-Mead
    polish, same convention as src/fso/optimization.optimize_instantaneous)."""
    eps_grid = np.linspace(0.01, 0.99, n_grid)
    mu_grid = np.linspace(0.01, 1.0, n_grid)
    EPS, MU = np.meshgrid(eps_grid, mu_grid, indexing="ij")
    R = sns_key_rate_eta(EPS, MU, eta_value, e_a, decoy_ratio=decoy_ratio, f_ec=f_ec, p_dark=p_dark)
    idx = np.unravel_index(np.argmax(R), R.shape)
    best_eps, best_mu, best_R = EPS[idx], MU[idx], R[idx]

    def neg_rate(x):
        e, m = x
        if not (0.01 <= e <= 0.99) or not (0.01 <= m <= 1.0):
            return 1e3
        return -float(sns_key_rate_eta(e, m, eta_value, e_a, decoy_ratio=decoy_ratio, f_ec=f_ec, p_dark=p_dark))

    res = minimize(neg_rate, x0=[best_eps, best_mu], method="Nelder-Mead",
                    options={"xatol": 1e-10, "fatol": 1e-14, "maxiter": 600})
    if res.success and -res.fun > best_R:
        return float(-res.fun)
    return float(best_R)


def R_deterministic(d_km, alpha_fso, eta_d, e_a, decoy_ratio, f_ec, p_dark):
    eta = float(fso_deterministic_eta(d_km, alpha=alpha_fso, eta_d=eta_d))
    return optimize_R_full(eta, e_a, decoy_ratio, f_ec, p_dark)


def R_gamma_gamma(d_km, alpha_fso, eta_d, e_a, Cn2, seed=0):
    sigma2 = float(rytov_variance(d_km, Cn2))
    if not bool(sigma_R2_within_validity(sigma2)):
        return None, sigma2  # out of GG model validity - not a meaningful comparison
    out = joint_optimize_ensemble(d_km, e_a, sigma2, alpha_fso=alpha_fso, eta_d=eta_d, seed=seed)
    return out["R"], sigma2


def run():
    rows = []

    R0_det = R_deterministic(**{k: BASE[k] for k in
                                 ("d_km", "alpha_fso", "eta_d", "e_a", "decoy_ratio", "f_ec", "p_dark")})
    R0_gg, sigma2_0 = R_gamma_gamma(BASE["d_km"], BASE["alpha_fso"], BASE["eta_d"], BASE["e_a"], BASE["Cn2"])
    print(f"Baseline (d={BASE['d_km']} km): R_det={R0_det:.6e}  "
          f"R_gg={R0_gg:.6e} (sigma_R2={sigma2_0:.3f}, valid={sigma_R2_within_validity(sigma2_0)})")

    # --- (A) deterministic-channel parameters ------------------------------
    det_params = ["d_km", "alpha_fso", "eta_d", "e_a", "decoy_ratio", "f_ec", "p_dark"]
    for pname in det_params:
        for pct in PERTURBATIONS:
            kwargs = {k: BASE[k] for k in det_params}
            kwargs[pname] = BASE[pname] * (1.0 + pct)
            R = R_deterministic(**kwargs)
            pct_change_R = (R - R0_det) / abs(R0_det) * 100.0
            elasticity = pct_change_R / (pct * 100.0)
            rows.append(dict(model="Deterministic FSO", parameter=pname,
                              baseline_value=BASE[pname], perturbed_value=kwargs[pname],
                              param_pct_change=pct * 100.0, R_baseline=R0_det, R_perturbed=R,
                              R_pct_change=pct_change_R, elasticity=elasticity, note=""))

    # --- (B) turbulence parameter Cn2 (Gamma-Gamma ensemble) ---------------
    for pct in PERTURBATIONS:
        Cn2_p = BASE["Cn2"] * (1.0 + pct)
        R, sigma2 = R_gamma_gamma(BASE["d_km"], BASE["alpha_fso"], BASE["eta_d"], BASE["e_a"], Cn2_p)
        if R is None:
            rows.append(dict(model="Gamma-Gamma FSO", parameter="Cn2",
                              baseline_value=BASE["Cn2"], perturbed_value=Cn2_p,
                              param_pct_change=pct * 100.0, R_baseline=R0_gg, R_perturbed=np.nan,
                              R_pct_change=np.nan, elasticity=np.nan,
                              note=f"sigma_R2={sigma2:.2f} exceeds validity ceiling "
                                   f"{SIGMA_R2_VALIDITY_MAX} - excluded"))
            continue
        pct_change_R = (R - R0_gg) / abs(R0_gg) * 100.0
        elasticity = pct_change_R / (pct * 100.0)
        rows.append(dict(model="Gamma-Gamma FSO", parameter="Cn2",
                          baseline_value=BASE["Cn2"], perturbed_value=Cn2_p,
                          param_pct_change=pct * 100.0, R_baseline=R0_gg, R_perturbed=R,
                          R_pct_change=pct_change_R, elasticity=elasticity,
                          note=f"sigma_R2={sigma2:.2f}"))

    df = pd.DataFrame(rows)

    # Rank parameters by |elasticity| at the +-10% perturbations
    ranking = (df[np.isclose(df["param_pct_change"].abs(), 10.0)]
               .groupby("parameter")["elasticity"].apply(lambda s: s.abs().mean())
               .sort_values(ascending=False))
    print("\nSensitivity ranking (mean |elasticity| at +-10%):")
    for pname, val in ranking.items():
        print(f"  {pname:>12s}: {val:8.3f}")

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV} ({len(df)} rows)")
    return df


if __name__ == "__main__":
    run()