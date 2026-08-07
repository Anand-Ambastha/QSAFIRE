"""Part 8 — Performance-Metrics Dataframe.

For every distance (and, for the Gamma-Gamma models, every turbulence case)
we jointly optimise (eps, mu') and tabulate: gain Qxx, QBER Exx, single-photon
yield s1, phase error e1^ph, Z-basis gain/QBER Sz/Ez, key rate R, and channel
loss (dB). For the Gamma-Gamma models, R is reported both from the fast
quadrature estimate and from an independent N_MC=1e5 Monte Carlo run at the
*same* optimised operating point, together with the Monte Carlo 95%
confidence interval.
"""

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from src.engine import full_metrics_eta
from src.fso.config import TURBULENCE_CASES, E_A, SECURE_THRESHOLD
from src.fso.channel_models import (
    fso_deterministic_eta, channel_loss_dB, rytov_variance,
    gamma_gamma_params, classify_turbulence_regime, validity_flag,
    sigma_R2_within_validity, SIGMA_R2_VALIDITY_MAX,
)
from src.fso.optimization import optimize_instantaneous, joint_optimize_ensemble
from src.fso.quadrature import gauss_laguerre_ensemble_metrics
from src.fso.monte_carlo import monte_carlo_ensemble_rate


def build_deterministic_table(distances_km, e_a=E_A):
    rows = []
    for d in tqdm(distances_km, desc="Deterministic FSO"):
        eta = float(fso_deterministic_eta(d))
        R_opt, eps_opt, mu_opt = optimize_instantaneous(eta, e_a)
        m = full_metrics_eta(eps_opt, mu_opt, eta, e_a)
        rows.append(dict(model="Deterministic FSO", d_km=d, L_total_km=2 * d,
                          eta=eta, loss_dB=float(channel_loss_dB(d)),
                          eps=eps_opt, mu=mu_opt, **{k: v for k, v in m.items() if k != "eta"}))
    return pd.DataFrame(rows)


def build_gamma_gamma_table(distances_km, turbulence_cases=TURBULENCE_CASES, e_a=E_A, n_mc=100_000, seed=0):
    rows = []
    for case_name, cn2 in turbulence_cases.items():
        for d in tqdm(distances_km, desc=f"Gamma-Gamma [{case_name}]"):
            sigma2 = float(rytov_variance(d, cn2))
            alpha_g, beta_g = gamma_gamma_params(sigma2)
            regime = str(classify_turbulence_regime(sigma2))
            flag = str(validity_flag(sigma2))
            opt = joint_optimize_ensemble(d, e_a, sigma2, seed=seed)
            eps_opt, mu_opt = opt["eps"], opt["mu"]
            m_quad = gauss_laguerre_ensemble_metrics(eps_opt, mu_opt, d, e_a, sigma2, n_points=40)
            mc = monte_carlo_ensemble_rate(eps_opt, mu_opt, d, e_a, sigma2, n_mc=n_mc, rng=np.random.default_rng(seed))
            rows.append(dict(
                model=f"Gamma-Gamma [{case_name}]", turbulence_case=case_name, d_km=d, L_total_km=2 * d,
                Cn2=cn2, sigma_R2=sigma2, alpha_g=alpha_g, beta_g=beta_g, regime=regime, validity_flag=flag,
                eta_det=float(fso_deterministic_eta(d)), loss_dB=float(channel_loss_dB(d)),
                eps=eps_opt, mu=mu_opt, R_quad=m_quad["R"], R_mc=mc["mean"], R_mc_std=mc["std"],
                R_mc_ci95_lo=mc["ci95_lo"], R_mc_ci95_hi=mc["ci95_hi"], Qxx=m_quad["Qxx"], Exx=m_quad["Exx"],
                s1=m_quad["s1"], e1ph=m_quad["e1ph"], Sz=m_quad["Sz"], Ez=m_quad["Ez"], best_method=opt["best_method"],
                within_validity=bool(sigma_R2_within_validity(sigma2)),
            ))
    df = pd.DataFrame(rows)
    n_bad = int((~df["within_validity"]).sum())
    if n_bad:
        print(f"[build_gamma_gamma_table] WARNING: {n_bad}/{len(df)} rows have "
              f"sigma_R2 > {SIGMA_R2_VALIDITY_MAX} (outside the Gamma-Gamma model's "
              f"validity range - see channel_models.SIGMA_R2_VALIDITY_MAX). At these "
              f"rows alpha_g/beta_g diverge unphysically and R_quad/R_mc will read as "
              f"numerically identical to the deterministic channel; do not use them to "
              f"draw GG-vs-deterministic comparisons. Filter on within_validity before "
              f"computing headline statistics.")
    return df


def maximum_performance_summary(df_det, df_gg):
    """NOTE: Gamma-Gamma rows with sigma_R2 beyond SIGMA_R2_VALIDITY_MAX are
    excluded before summarising - at those distances alpha_g/beta_g diverge
    unphysically and the GG channel numerically collapses onto the
    deterministic one (see channel_models.gamma_gamma_params docstring), so
    including them would silently bias Max_SKR / Max_Secure_Distance toward
    "GG == deterministic" as an artifact of extrapolating the model past its
    valid range rather than a real result."""
    rows = []

    def _summarize(sub, model_name, r_col, loss_col="loss_dB", d_col="d_km", truncated=False):
        sub = sub.sort_values(d_col)
        i_max = sub[r_col].idxmax()
        max_R, d_at_max = sub.loc[i_max, r_col], sub.loc[i_max, d_col]
        secure = sub[sub[r_col] > SECURE_THRESHOLD]
        d_range_max = sub[d_col].max()
        if len(secure) > 0:
            max_secure_d = secure[d_col].max()
            max_loss = secure.loc[secure[d_col].idxmax(), loss_col]
            if max_secure_d >= d_range_max and sub.loc[sub[d_col].idxmax(), r_col] > SECURE_THRESHOLD:
                note = "secure throughout simulated range - true cutoff lies beyond it"
            else:
                note = ""
        else:
            max_secure_d, max_loss, note = np.nan, np.nan, "threshold not crossed anywhere in simulated distance range"
        if truncated:
            note = (note + "; " if note else "") + \
                f"distance range truncated to sigma_R2<={SIGMA_R2_VALIDITY_MAX} (GG model validity)"
        rows.append(dict(Model=model_name, Max_SKR=max_R, Distance_at_Max_SKR_km=d_at_max,
                          Max_Secure_Distance_km=max_secure_d, Max_Tolerable_Loss_dB=max_loss, Note=note))

    _summarize(df_det, "Deterministic FSO", "R")
    for case in df_gg["turbulence_case"].unique():
        sub_full = df_gg[df_gg["turbulence_case"] == case]
        sub_valid = sub_full[sub_full["within_validity"]] if "within_validity" in sub_full.columns else sub_full
        was_truncated = len(sub_valid) < len(sub_full)
        if sub_valid.empty:
            rows.append(dict(Model=f"Gamma-Gamma MC [{case}]", Max_SKR=np.nan, Distance_at_Max_SKR_km=np.nan,
                              Max_Secure_Distance_km=np.nan, Max_Tolerable_Loss_dB=np.nan,
                              Note=f"no rows within sigma_R2<={SIGMA_R2_VALIDITY_MAX} validity range"))
            continue
        _summarize(sub_valid, f"Gamma-Gamma MC [{case}]", "R_mc", truncated=was_truncated)
        _summarize(sub_valid, f"Gamma-Gamma Quadrature [{case}]", "R_quad", truncated=was_truncated)
    return pd.DataFrame(rows)