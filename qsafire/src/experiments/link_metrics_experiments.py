"""Part XIII orchestrator: Per-Link SNS-TF-QKD Metrics (Phase 9 final).

Reproduces notebook cells 238-249 in order. Consumes ``links_state`` (Part
XII: ``pass_results``, ``link_eval_cache``). Returns a ``link_metrics_state``
dict.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from src.pat.losses import E_A_MISALIGNMENT, R_TX_HZ, R_NOISE_BASELINE_HZ, P_DARK_REAL
from src.links.city_pairs import LINK_PAIRS
from src.links.metrics import (
    CITY_BG_MULTIPLIER, validate_symmetric_limit, optimize_joint_asym_ensemble, sns_key_rate_asym_full,
)

MARGIN_FACTOR = 3.0


def per_link_breakeven_and_metrics(pass_results, link_eval_cache):
    """Cell 244."""
    link_metrics_asym = {}
    for city_a, city_b in LINK_PAIRS:
        ev = link_eval_cache[(city_a, city_b)]
        eta_A, eta_B = ev['A']['eta_total'], ev['B']['eta_total']
        sigma_A, sigma_B = ev['A']['sigma_R2'], ev['B']['sigma_R2']

        p_bg_A = CITY_BG_MULTIPLIER[city_a]['mult'] * R_NOISE_BASELINE_HZ / R_TX_HZ
        p_bg_B = CITY_BG_MULTIPLIER[city_b]['mult'] * R_NOISE_BASELINE_HZ / R_TX_HZ
        p_bg_combined_ref = p_bg_A + p_bg_B

        def R_of_supp(supp, eta_A=eta_A, eta_B=eta_B, sigma_A=sigma_A, sigma_B=sigma_B, p_bg=p_bg_combined_ref):
            p_dark_test = P_DARK_REAL + p_bg/supp
            R, _, _, _ = optimize_joint_asym_ensemble(eta_A, eta_B, sigma_A, sigma_B, E_A_MISALIGNMENT,
                                                        p_dark_test, n_mc=3000, n_grid=5)
            return R

        try:
            thresh = brentq(R_of_supp, 1.0, 5000.0, xtol=2.0)
        except ValueError:
            thresh = np.nan
        adopted = MARGIN_FACTOR*thresh if np.isfinite(thresh) else np.nan
        p_dark_final = P_DARK_REAL + p_bg_combined_ref/adopted if np.isfinite(adopted) else P_DARK_REAL + p_bg_combined_ref

        R_f, eps_f, mu2_f, mu1_f = optimize_joint_asym_ensemble(eta_A, eta_B, sigma_A, sigma_B, E_A_MISALIGNMENT,
                                                                  p_dark_final, n_mc=10000, n_grid=8)
        full = sns_key_rate_asym_full(eps_f, mu2_f, mu1_f, eta_A, eta_B, E_A_MISALIGNMENT, p_dark=p_dark_final)
        full = {k: float(np.asarray(v)) for k, v in full.items()}

        link_metrics_asym[(city_a, city_b)] = dict(
            eta_A=eta_A, eta_B=eta_B, sigma_A=sigma_A, sigma_B=sigma_B,
            breakeven_suppression=thresh, adopted_suppression=adopted, p_dark_final=p_dark_final,
            eps=eps_f, mu2_signal=mu2_f, mu1_decoy=mu1_f, **full)

        print(f"{city_a:10s}-{city_b:10s}: eta_A={eta_A:.3e} eta_B={eta_B:.3e}  breakeven={thresh:7.1f}x  "
              f"adopted={adopted:7.1f}x  QBER={full['Ez']*100:.2f}%  phase_err={full['e1ph']*100:.2f}%  "
              f"R={full['R']:.3e} bits/pulse ({max(full['R'], 0)*R_TX_HZ/1000:.4f} kbit/s)")

    return link_metrics_asym


def master_link_table(pass_results, link_metrics_asym):
    """Cell 246."""
    rows = []
    for (city_a, city_b), m in link_metrics_asym.items():
        rows.append(dict(
            link=f"{city_a}-{city_b}", distance_km=round(pass_results[(city_a, city_b)]['dist_km'], 1),
            eta_A=f"{m['eta_A']:.2e}", eta_B=f"{m['eta_B']:.2e}",
            sigma_A=round(m['sigma_A'], 3), sigma_B=round(m['sigma_B'], 3),
            breakeven_supp=round(m['breakeven_suppression'], 1) if np.isfinite(m['breakeven_suppression']) else "n/a",
            adopted_supp=round(m['adopted_suppression'], 1) if np.isfinite(m['adopted_suppression']) else "n/a",
            eps_star=round(m['eps'], 4), mu2_signal=round(m['mu2_signal'], 4), mu1_decoy=round(m['mu1_decoy'], 5),
            gain_Sz=f"{m['Sz']:.3e}", QBER_Ez_pct=round(m['Ez']*100, 2), phase_err_e1ph_pct=round(m['e1ph']*100, 2),
            yield_s1=f"{m['s1']:.3e}", SKR_bits_per_pulse=f"{m['R']:.3e}",
            SKR_kbps=round(max(m['R'], 0)*R_TX_HZ/1000, 4),
        ))
    df_asym = pd.DataFrame(rows)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 230)
    print(df_asym.to_string(index=False))
    return df_asym


def plot_final_metrics_summary(link_metrics_asym, save_path="phase13final_metrics_summary.png", show=True):
    """Cell 247."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    labels = [f"{a}-{b}" for a, b in LINK_PAIRS]
    bar_colors = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]

    axes[0].bar(labels, [link_metrics_asym[p]['Ez']*100 for p in LINK_PAIRS], color=bar_colors)
    axes[0].set_ylabel('QBER Ez (%)'); axes[0].set_title('QBER by link (rigorous asymmetric)')
    axes[0].tick_params(axis='x', rotation=20)

    axes[1].bar(labels, [link_metrics_asym[p]['e1ph']*100 for p in LINK_PAIRS], color=bar_colors)
    axes[1].set_ylabel('Phase-error rate e1_ph (%)'); axes[1].set_title('Phase error by link')
    axes[1].tick_params(axis='x', rotation=20)

    axes[2].bar(labels, [max(link_metrics_asym[p]['R'], 0)*R_TX_HZ/1000 for p in LINK_PAIRS], color=bar_colors)
    axes[2].set_ylabel('SKR (kbit/s)'); axes[2].set_title('Secure key rate by link')
    axes[2].tick_params(axis='x', rotation=20)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


def run_full_part13_pipeline(links_state, show_plots=True):
    """Runs notebook cells 238-249 in order. Returns link_metrics_state dict."""
    validate_symmetric_limit()
    link_metrics_asym = per_link_breakeven_and_metrics(links_state['pass_results'], links_state['link_eval_cache'])
    df_asym = master_link_table(links_state['pass_results'], link_metrics_asym)
    plot_final_metrics_summary(link_metrics_asym, show=show_plots)
    return dict(link_metrics_asym=link_metrics_asym, df_asym=df_asym)
