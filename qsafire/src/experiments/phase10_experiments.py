"""Part X orchestrator: Joint Signal/Decoy Intensity and Sending-Probability
Optimization (Phase 8).

Reproduces notebook cells 196-212 in order. Consumes ``state`` (Part III),
``phase5_state`` (Part VII), ``turb_state`` (Part V), ``phase6_state``
(Part VIII, including REFERENCE_BITS_PER_PASS-equivalent data and
p_dark_suppressed). Returns a ``phase10_state`` dict.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from src.satellite.geometry import SITES
from src.pat.losses import E_A_MISALIGNMENT, R_TX_HZ, P_DARK_REAL
from src.coupling.phase6 import optimize_ensemble_real
from src.experiments.phase6_experiments import SUBSAMPLE_STRIDE
from src.optimization.joint import (
    sns_key_rate_eta_joint, refactor_correctness_check, optimize_joint_signal_decoy,
    optimize_joint_ensemble,
)

COLORS = {"Delhi (Alice)": "#1f77b4", "Mumbai (Bob)": "#d62728"}


def quantify_improvement(state, phase5_results, turb_results_site, p_dark_suppressed):
    """Cell 204: joint vs. fixed-decoy optimizer, at each station's
    best-elevation point."""
    idx = state['idx']
    improvement_rows = []
    for station in SITES:
        eta_best = phase5_results[station]['eta_total'][idx].max()
        sigma_best = float(turb_results_site[station]['sigma_R2'][idx].max())
        p_dark_val = p_dark_suppressed[station]

        R_old, eps_old, mu_old = optimize_ensemble_real(eta_best, E_A_MISALIGNMENT, sigma_best, p_dark_val,
                                                          n_mc=20000, n_grid=12)
        R_new, eps_new, mu2_new, mu1_new = optimize_joint_ensemble(eta_best, E_A_MISALIGNMENT, sigma_best,
                                                                     p_dark_val, n_mc=20000, n_grid=10)
        pct = 100*(R_new-R_old)/abs(R_old) if R_old != 0 else np.nan
        improvement_rows.append(dict(station=station, R_old_fixed_decoy=R_old, eps_old=eps_old, mu_old=mu_old,
                                      R_new_joint=R_new, eps_new=eps_new, mu2_new=mu2_new, mu1_new=mu1_new,
                                      decoy_ratio_found=mu1_new/mu2_new if mu2_new else np.nan,
                                      pct_improvement=pct))
    return pd.DataFrame(improvement_rows)


def rerun_pass_with_joint_optimizer(state, phase5_results, turb_results_site, p_dark_suppressed,
                                     REFERENCE_BITS_PER_PASS):
    """Cell 207: repeat Phase 6 Section 8.6's per-time-slice integration
    exactly, with optimize_joint_ensemble in place of optimize_ensemble_real."""
    t_grid, idx = state['t_grid'], state['idx']
    dt = t_grid[1] - t_grid[0]
    idx_sub = idx[::SUBSAMPLE_STRIDE]
    if idx_sub[-1] != idx[-1]:
        idx_sub = np.append(idx_sub, idx[-1])
    t_sub = t_grid[idx_sub]

    phase10_results = {}
    for station in SITES:
        R_sub = np.full(len(idx_sub), np.nan)
        for k, i in enumerate(idx_sub):
            eta_i = phase5_results[station]['eta_total'][i]
            sigma_R2_i = float(turb_results_site[station]['sigma_R2'][i])
            R_opt, _, _, _ = optimize_joint_ensemble(eta_i, E_A_MISALIGNMENT, sigma_R2_i,
                                                       p_dark_suppressed[station], n_mc=8000, n_grid=8)
            R_sub[k] = R_opt
        R_full = np.interp(t_grid[idx], t_sub, R_sub)
        phase10_results[station] = dict(t=t_grid[idx], R=R_full)

    print("Total key bits per pass, OLD (fixed decoy ratio) vs. NEW (joint signal+decoy optimization):\n")
    for station in SITES:
        old_bits = REFERENCE_BITS_PER_PASS[station]
        new_bits = np.sum(np.clip(phase10_results[station]['R'], 0, None)) * R_TX_HZ * dt
        print(f"  {station:16s}: OLD = {old_bits:.3e} bits/pass ({old_bits/(len(idx)*dt)/1000:.3f} kbit/s avg)   "
              f"NEW = {new_bits:.3e} bits/pass ({new_bits/(len(idx)*dt)/1000:.3f} kbit/s avg)   "
              f"({100*(new_bits-old_bits)/old_bits:+.1f}%)")

    return phase10_results


def plot_joint_opt_comparison(phase6_results, phase10_results, save_path="phase10_joint_opt_comparison.png", show=True):
    """Cell 208."""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for station in SITES:
        ax.plot(phase6_results[station]['t'], np.clip(phase6_results[station]['R'], 0, None)*R_TX_HZ/1000,
                color=COLORS[station], lw=1.3, ls='--', label=f'{station} (old, fixed decoy)')
        ax.plot(phase10_results[station]['t'], np.clip(phase10_results[station]['R'], 0, None)*R_TX_HZ/1000,
                color=COLORS[station], lw=2.0, label=f'{station} (new, joint opt.)')
    ax.set_xlabel('Time relative to engineered crossing (s)')
    ax.set_ylabel('Key rate (kbit/s)')
    ax.set_title('Effect of joint signal+decoy optimization on the per-pass key rate')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


def recheck_suppression_breakeven(state, phase5_results, turb_results_site, suppression_needed):
    """Cell 210: recompute the breakeven suppression factor with the joint
    optimizer in place of optimize_ensemble_real."""
    idx = state['idx']
    new_suppression_needed = {}
    for station in SITES:
        eta_best = phase5_results[station]['eta_total'][idx].max()
        p_bg_current = phase5_results[station]['p_dark_total'] - P_DARK_REAL
        sigma_R2_here = float(turb_results_site[station]['sigma_R2'][idx].max())

        def R_of_suppression(supp, p_bg_current=p_bg_current, eta_best=eta_best, sigma_R2_here=sigma_R2_here):
            p_dark_test = P_DARK_REAL + p_bg_current/supp
            R, _, _, _ = optimize_joint_ensemble(eta_best, E_A_MISALIGNMENT, sigma_R2_here, p_dark_test,
                                                   n_mc=8000, n_grid=8)
            return R

        try:
            supp_thresh = brentq(R_of_suppression, 1.0, 450.0, xtol=0.5)
        except ValueError:
            supp_thresh = float('nan')
        new_suppression_needed[station] = supp_thresh
        old_val = suppression_needed[station]
        print(f"{station:16s}: NEW breakeven suppression = {supp_thresh:.1f}x   "
              f"(OLD, fixed-decoy requirement was {old_val:.0f}x -- "
              f"{100*(supp_thresh-old_val)/old_val:+.1f}%)")

    return new_suppression_needed


def run_full_part10_pipeline(state, phase5_state, turb_state, phase6_state, weather_state, show_plots=True):
    """Runs notebook cells 196-212 in order. Returns phase10_state dict."""
    refactor_correctness_check()

    phase5_results = phase5_state['phase5_results']
    turb_results_site = turb_state['turb_results_site']
    p_dark_suppressed = phase6_state['p_dark_suppressed']
    suppression_needed = phase6_state['suppression_needed']
    phase6_results = phase6_state['phase6_results']
    REFERENCE_BITS_PER_PASS = weather_state['REFERENCE_BITS_PER_PASS']

    df_improve = quantify_improvement(state, phase5_results, turb_results_site, p_dark_suppressed)
    print(df_improve)

    phase10_results = rerun_pass_with_joint_optimizer(state, phase5_results, turb_results_site,
                                                        p_dark_suppressed, REFERENCE_BITS_PER_PASS)
    plot_joint_opt_comparison(phase6_results, phase10_results, show=show_plots)
    new_suppression_needed = recheck_suppression_breakeven(state, phase5_results, turb_results_site,
                                                             suppression_needed)

    return dict(df_improve=df_improve, phase10_results=phase10_results,
                new_suppression_needed=new_suppression_needed)
