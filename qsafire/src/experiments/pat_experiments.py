"""Part VII orchestrator: PAT, Optics and Detector Real-Loss Module (Phase 5).

Reproduces notebook cells 143-161 in order. Consumes ``state`` (Part III) and
``atm_results`` (Part IV baseline). Returns a ``phase5_state`` dict carrying
``phase5_results`` forward to Part VIII.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.satellite.geometry import SITES
from src.pat.losses import (
    LAMBDA_P5_M, THETA_DIV_RAD, D_RX_M, OBSCURATION_RATIO, D_RX_INT_M, R_TX_HZ,
    G_TX, G_RX, ETA_RX_INTERNAL, POINTING_CASES, ETA_DET_REAL, DARK_COUNT_RATE_HZ,
    TIMING_JITTER_NS, AFTERPULSE_PROB, P_DARK_REAL, R_NOISE_BASELINE_HZ,
    BACKGROUND_SCENARIOS, E_A_MISALIGNMENT, eta_diffraction, eta_pointing,
    optimize_instantaneous_real,
)

COLORS = {"Delhi (Alice)": "#1f77b4", "Mumbai (Bob)": "#d62728"}


def validate_diffraction():
    """Cell 145 (script part)."""
    L_test = 1200e3
    loss_test_dB = -10*np.log10(eta_diffraction(L_test))
    print(f"G_Tx = {G_TX:.3e}, G_Rx = {G_RX:.3e}")
    print(f"Computed diffraction loss @ 1200 km: {loss_test_dB:.2f} dB  "
          f"(Liao et al. 2017 report ~22 dB at this range -- same order, aperture/divergence details differ slightly)")
    print(f"Receiver internal loss (7 mirrors, 95% each): {10*np.log10(ETA_RX_INTERNAL):.2f} dB")


def pointing_table():
    """Cell 147 (script part)."""
    for label, p in POINTING_CASES.items():
        eta_p = eta_pointing(p['sigma_j'])
        print(f"{label:32s}: sigma_jitter={p['sigma_j']*1e6:.2f} urad -> eta_pointing={eta_p:.4f} "
              f"({10*np.log10(eta_p):+.3f} dB)   [{p['conf']}]")


def detector_realism_summary():
    """Cell 149 (script part)."""
    print(f"Real per-pulse dark-count probability: {P_DARK_REAL:.2e}  "
          f"(vs. the notebook's idealized default P_DARK=1e-11 -- {P_DARK_REAL/1e-11:.0f}x higher / more realistic)")
    print(f"Real detector efficiency: {ETA_DET_REAL} (vs. idealized default ETA_D=0.8)")
    print(f"Timing jitter < {TIMING_JITTER_NS} ns, afterpulsing < {AFTERPULSE_PROB*100:.2f}% "
          f"-- not separately modelled in the key-rate formula (Wang/Yu/Hu Eq. (4) does not have timing-"
          f"jitter or afterpulsing terms); carried here as reported hardware characteristics only.")


def background_table():
    """Cell 151 (script part)."""
    rows = []
    for label, p in BACKGROUND_SCENARIOS.items():
        R_noise = R_NOISE_BASELINE_HZ * p['mult']
        p_bg = R_noise / R_TX_HZ
        rows.append(dict(scenario=label, R_noise_Hz=R_noise, p_background_per_pulse=p_bg,
                          p_dark_plus_background=p_bg + P_DARK_REAL, confidence=p['conf']))
    return pd.DataFrame(rows)


def run_real_loss_pass(state, atm_results):
    """Cell 153 (script part): build the fully-dressed real channel and
    re-run the SNS-TF-QKD engine with real losses, night scenario."""
    results, idx = state['results'], state['idx']
    sigma_j_case = POINTING_CASES["Conservative default ([S2025])"]['sigma_j']
    eta_point_val = eta_pointing(sigma_j_case)

    phase5_results = {}
    for station in SITES:
        L_m = results[station]['rng'] * 1000.0
        eta_diff_t = eta_diffraction(L_m)
        eta_atm_t = atm_results[station]['eta_atm']            # Phase 3 baseline
        eta_total_t = eta_diff_t * eta_atm_t * eta_point_val * ETA_RX_INTERNAL * ETA_DET_REAL

        bg_key = f"{station.split(' ')[0]} ({'Alice' if 'Delhi' in station else 'Bob'}) -- night"
        p_bg = BACKGROUND_SCENARIOS[bg_key]['mult'] * R_NOISE_BASELINE_HZ / R_TX_HZ
        p_dark_total = P_DARK_REAL + p_bg

        R_t = np.full_like(eta_total_t, np.nan)
        eps_t = np.full_like(eta_total_t, np.nan)
        mu_t = np.full_like(eta_total_t, np.nan)
        for i in idx:  # only optimize within the common dual-visibility window (rest is not usable anyway)
            R_opt, eps_opt, mu_opt = optimize_instantaneous_real(eta_total_t[i], E_A_MISALIGNMENT, p_dark_total)
            R_t[i], eps_t[i], mu_t[i] = R_opt, eps_opt, mu_opt

        phase5_results[station] = dict(eta_diff=eta_diff_t, eta_total=eta_total_t,
                                        p_dark_total=p_dark_total, R=R_t, eps=eps_t, mu=mu_t)

    print("Real-loss key rate (bits/pulse) within the common dual-visibility window, night scenario:")
    for station in SITES:
        R_win = phase5_results[station]['R'][idx]
        valid = R_win[~np.isnan(R_win)]
        positive = valid[valid > 0]
        print(f"  {station:16s}: max R = {valid.max():.3e} bits/pulse   "
              f"positive-key fraction of window = {100*len(positive)/len(valid):.0f}%   "
              f"p_dark_total = {phase5_results[station]['p_dark_total']:.2e}")

    return phase5_results


def root_cause_diagnostic(phase5_results, idx):
    """Cell 155."""
    print("Isolating the cause, at each station's best (lowest-loss) point in the window, e_a=0:\n")
    for station in SITES:
        eta_best = phase5_results[station]['eta_total'][idx].max()
        p_dark_val = phase5_results[station]['p_dark_total']
        R_asis, _, _ = optimize_instantaneous_real(eta_best, 0.0, p_dark_val)
        R_dark_only, _, _ = optimize_instantaneous_real(eta_best, 0.0, P_DARK_REAL)  # strip out background entirely
        print(f"{station}: best-point loss = {-10*np.log10(eta_best):.2f} dB")
        print(f"  as computed (dark + background, e_a=0):        R = {R_asis:.3e}  ({'secure' if R_asis>0 else 'NO KEY'})")
        print(f"  dark-count only, background removed:            R = {R_dark_only:.3e}  ({'secure' if R_dark_only>0 else 'NO KEY'})")
        print(f"  -> background light alone flips this station from a viable link to no link.\n")


def plot_real_loss_keyrate(state, phase5_results, save_path="phase5_real_loss_keyrate.png", show=True):
    """Cell 157."""
    t_grid, t_start, t_end = state['t_grid'], state['t_start'], state['t_end']
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for station in SITES:
        axes[0].plot(t_grid, -10*np.log10(np.clip(phase5_results[station]['eta_total'], 1e-30, 1)),
                     color=COLORS[station], lw=1.8, label=station)
    axes[0].axvspan(t_start, t_end, color='green', alpha=0.12, label='common dual-visibility window')
    axes[0].set_ylabel('Total real-loss channel loss (dB)')
    axes[0].set_title('Fully-dressed channel loss (diffraction + atmosphere + pointing + optics + detector)')
    axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3)

    for station in SITES:
        R_plot = phase5_results[station]['R'].copy()
        axes[1].plot(t_grid, R_plot, color=COLORS[station], lw=1.8, label=station)
    axes[1].axhline(0, color='gray', lw=0.8)
    axes[1].axvspan(t_start, t_end, color='green', alpha=0.12)
    axes[1].set_ylabel('Key rate R (bits/pulse)')
    axes[1].set_xlabel('Time relative to engineered crossing (s)')
    axes[1].set_title('Real-loss SNS-TF-QKD key rate during the pass (night scenario)')
    axes[1].legend(fontsize=9); axes[1].grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


def day_vs_night(phase5_results, idx):
    """Cell 159."""
    for station in SITES:
        eta_best = phase5_results[station]['eta_total'][idx].max()
        p_bg_day = BACKGROUND_SCENARIOS["Either city -- day"]['mult'] * R_NOISE_BASELINE_HZ / R_TX_HZ
        p_dark_day = P_DARK_REAL + p_bg_day

        R_dark_only, _, _ = optimize_instantaneous_real(eta_best, 0.0, P_DARK_REAL)
        R_night, _, _ = optimize_instantaneous_real(eta_best, 0.0, phase5_results[station]['p_dark_total'])
        R_day, _, _ = optimize_instantaneous_real(eta_best, 0.0, p_dark_day)
        print(f"{station} @ best-elevation point in the pass (e_a=0):")
        print(f"  dark-count only (no background):  R = {R_dark_only:.3e}  ({'secure' if R_dark_only > 0 else 'NO KEY'})")
        print(f"  + night urban background:         R = {R_night:.3e}  ({'secure' if R_night > 0 else 'NO KEY'})")
        print(f"  + day background (x{BACKGROUND_SCENARIOS['Either city -- day']['mult']:.0f} vs. night baseline): "
              f"R = {R_day:.3e}  ({'secure' if R_day > 0 else 'NO KEY'})\n")


def run_full_part7_pipeline(state, atm_results, show_plots=True):
    """Runs notebook cells 143-161 in order. Returns phase5_state dict for
    Part VIII onward."""
    validate_diffraction()
    pointing_table()
    detector_realism_summary()
    df_background = background_table()
    print(df_background)
    phase5_results = run_real_loss_pass(state, atm_results)
    root_cause_diagnostic(phase5_results, state['idx'])
    plot_real_loss_keyrate(state, phase5_results, show=show_plots)
    day_vs_night(phase5_results, state['idx'])
    return dict(phase5_results=phase5_results)
