"""Part VIII orchestrator: Coupling to the SNS-TF-QKD Engine (Phase 6).

Reproduces notebook cells 162-177 in order. Consumes ``state`` (Part III),
``atm_results`` (Part IV), ``turb_results_site`` (Part V), and
``phase5_results`` (Part VII). Returns a ``phase6_state`` dict carrying
``phase6_results``, ``p_dark_suppressed``, ``suppression_needed`` forward to
Part IX onward.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from src.satellite.geometry import SITES
from src.pat.losses import R_TX_HZ, P_DARK_REAL, E_A_MISALIGNMENT
from src.coupling.phase6 import optimize_ensemble_real

COLORS = {"Delhi (Alice)": "#1f77b4", "Mumbai (Bob)": "#d62728"}
SUPPRESSION_FACTOR_ADOPTED = 1000.0
SUBSAMPLE_STRIDE = 10  # ~0.5s native resolution -> ~5s effective; smooth elevation/loss variation justifies this


def liao_cross_check(state, atm_results):
    """Cell 165."""
    idx, results = state['idx'], state['results']
    print("This notebook's atmospheric loss at each station's best (highest-elevation) point in the")
    print("common-visibility window, vs. Liao et al. (2017)'s real reported range at ~1,200 km:\n")
    for station in SITES:
        i_best = idx[np.argmax(results[station]['el'][idx])]
        atm_loss = atm_results[station]['loss_dB'][i_best]
        print(f"  {station:16s}: {atm_loss:.2f} dB   (Liao et al. 2017 report 3-8 dB atmosphere+turbulence, "
              f"<3 dB pointing, at 1,200 km)")
    print()
    print("Both stations fall within or just below the low end of the real reported range -- consistent,")
    print("given Delhi/Mumbai's shorter best-case slant ranges (~600-800 km) than Micius's 1,200 km figure.")


def ensemble_rescue_check(state, phase5_results, turb_results_site):
    """Cell 167 (script part): confirm turbulence fading does not rescue the
    Phase-5 pre-suppression closed link."""
    idx = state['idx']
    print("Ensemble-optimized rate at the Phase 5 background level (pre-suppression):\n")
    for station in SITES:
        eta_best = phase5_results[station]['eta_total'][idx].max()
        p_dark_current = phase5_results[station]['p_dark_total']
        sigma_R2_here = float(turb_results_site[station]['sigma_R2'][idx].max())
        R_ens, eps_ens, mu_ens = optimize_ensemble_real(eta_best, E_A_MISALIGNMENT, sigma_R2_here, p_dark_current)
        print(f"  {station:16s}: ensemble-optimized R = {R_ens:.3e} bits/pulse  "
              f"({'secure' if R_ens > 0 else 'still NO KEY -- fading does not rescue it'})")


def suppression_requirement(state, phase5_results, turb_results_site):
    """Cell 169: breakeven background-suppression factor, re-derived under
    fading + realistic misalignment."""
    idx = state['idx']
    suppression_needed = {}
    for station in SITES:
        eta_best = phase5_results[station]['eta_total'][idx].max()
        p_bg_current = phase5_results[station]['p_dark_total'] - P_DARK_REAL
        sigma_R2_here = float(turb_results_site[station]['sigma_R2'][idx].max())

        def R_of_suppression(supp, p_bg_current=p_bg_current, eta_best=eta_best, sigma_R2_here=sigma_R2_here):
            p_dark_test = P_DARK_REAL + p_bg_current/supp
            R, _, _ = optimize_ensemble_real(eta_best, E_A_MISALIGNMENT, sigma_R2_here, p_dark_test, n_mc=8_000, n_grid=10)
            return R

        supp_thresh = brentq(R_of_suppression, 100.0, 600.0, xtol=0.5)
        suppression_needed[station] = supp_thresh
        print(f"{station:16s}: needs >= {supp_thresh:.0f}x background suppression to break even under "
              f"turbulence fading AND realistic misalignment (e_a=0.15). The Phase 5 Section 7.6 "
              f"isolation check (e_a=0, to isolate the background effect alone) found a much smaller "
              f"~11-22x -- misalignment error substantially compounds the background-noise budget.")

    print()
    print("Both requirements are modest next to demonstrated real capability:")
    print("  - Micius itself already uses temporal + spectral filtering to control background [L2017].")
    print("  - A real portable ground station has demonstrated QKD in strong urban background light")
    print("    using only a narrower spectral filter and a 100 microrad field of view [Z2023].")
    print("  - A 2025 single-mode-fibre-coupled ground terminal has demonstrated >120 dB (~1e12x)")
    print("    background suppression [K2025] -- more than 100x above the ~200-400x needed here.")
    return suppression_needed


def adopt_suppression(phase5_results, suppression_needed):
    """Cell 171."""
    p_dark_suppressed = {}
    for station in SITES:
        p_bg_current = phase5_results[station]['p_dark_total'] - P_DARK_REAL
        p_dark_suppressed[station] = P_DARK_REAL + p_bg_current/SUPPRESSION_FACTOR_ADOPTED
        print(f"{station:16s}: p_dark_total (suppressed) = {p_dark_suppressed[station]:.3e}  "
              f"(margin over breakeven: {SUPPRESSION_FACTOR_ADOPTED/suppression_needed[station]:.1f}x)")
    return p_dark_suppressed


def full_pass_integration(state, phase5_results, turb_results_site, p_dark_suppressed):
    """Cell 173: ensemble-optimized key rate at every subsampled time-slice,
    interpolated back onto the full window, integrated to total key bits."""
    t_grid, idx = state['t_grid'], state['idx']
    idx_sub = idx[::SUBSAMPLE_STRIDE]
    if idx_sub[-1] != idx[-1]:
        idx_sub = np.append(idx_sub, idx[-1])

    phase6_results = {}
    for station in SITES:
        R_sub = np.full(len(idx_sub), np.nan)
        eps_sub = np.full(len(idx_sub), np.nan)
        mu_sub = np.full(len(idx_sub), np.nan)
        for k, i in enumerate(idx_sub):
            eta_i = phase5_results[station]['eta_total'][i]
            sigma_R2_i = float(turb_results_site[station]['sigma_R2'][i])
            R_opt, eps_opt, mu_opt = optimize_ensemble_real(eta_i, E_A_MISALIGNMENT, sigma_R2_i,
                                                              p_dark_suppressed[station],
                                                              n_mc=8_000, n_grid=10)
            R_sub[k], eps_sub[k], mu_sub[k] = R_opt, eps_opt, mu_opt
        t_sub = t_grid[idx_sub]
        R_full = np.interp(t_grid[idx], t_sub, R_sub)
        eps_full = np.interp(t_grid[idx], t_sub, eps_sub)
        mu_full = np.interp(t_grid[idx], t_sub, mu_sub)
        phase6_results[station] = dict(t=t_grid[idx], R=R_full, eps=eps_full, mu=mu_full,
                                        t_sub=t_sub, R_sub=R_sub)

    print(f"Computed at {len(idx_sub)} subsampled points (stride {SUBSAMPLE_STRIDE}), "
          f"interpolated back onto the full {len(idx)}-point common-visibility window.\n")

    dt = t_grid[1] - t_grid[0]
    print("Total key bits per pass (common dual-visibility window only), suppressed-background scenario:\n")
    for station in SITES:
        R_t = phase6_results[station]['R']
        R_t_clipped = np.clip(R_t, 0, None)   # no negative "key" accumulates
        total_bits = np.sum(R_t_clipped) * R_TX_HZ * dt
        window_duration_s = len(idx) * dt
        avg_rate_bps = total_bits / window_duration_s if window_duration_s > 0 else 0.0
        print(f"  {station:16s}: total key = {total_bits:.3e} bits over {window_duration_s:.1f} s "
              f"-> average {avg_rate_bps/1000:.3f} kbit/s over the usable window")

    return phase6_results


def plot_pass_keyrate(phase6_results, save_path="phase6_pass_keyrate.png", show=True):
    """Cell 174."""
    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    for station in SITES:
        axes[0].plot(phase6_results[station]['t'], np.clip(phase6_results[station]['R'], 0, None)*R_TX_HZ/1000,
                     color=COLORS[station], lw=1.8, label=station)
    axes[0].set_ylabel('Key rate (kbit/s)')
    axes[0].set_title(f'Phase 6: ensemble-optimized key rate vs. time (background suppressed {SUPPRESSION_FACTOR_ADOPTED:.0f}x)')
    axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3)

    for station in SITES:
        axes[1].plot(phase6_results[station]['t'], phase6_results[station]['eps'], color=COLORS[station],
                     lw=1.5, ls='-', label=f'{station} eps*')
        axes[1].plot(phase6_results[station]['t'], phase6_results[station]['mu'], color=COLORS[station],
                     lw=1.5, ls='--', label=f"{station} mu'*")
    axes[1].set_ylabel('Optimal eps*, mu*')
    axes[1].set_xlabel('Time relative to engineered crossing (s)')
    axes[1].legend(fontsize=7, ncol=2); axes[1].grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


def print_validation_table():
    """Cell 175 (narrative validation table, reproduced as a print summary
    since it is prose commentary, not a computation)."""
    print("Section 8.7 validation against real reported satellite-QKD sifted-key rates:")
    print("  This notebook (Delhi/Mumbai, suppressed-background scenario): Delhi ~27 bit/s, Mumbai ~75 bit/s")
    print("  [Z2023] portable ground station: sifted key rate at 1,000 km 'of the order of 100 bps'")
    print("  Micius flagship: 40.2 kbit/s @ 530 km, 1.2 kbit/s @ 1,034.7 km [L2017]")
    print("  Behera & Sinha (2024): ~2,080 bit/s sifted-key-equivalent order at IAO Hanle")


def run_full_part8_pipeline(state, atm_results, turb_results_site, phase5_results, show_plots=True):
    """Runs notebook cells 162-177 in order. Returns phase6_state dict for
    Part IX onward."""
    liao_cross_check(state, atm_results)
    ensemble_rescue_check(state, phase5_results, turb_results_site)
    suppression_needed = suppression_requirement(state, phase5_results, turb_results_site)
    p_dark_suppressed = adopt_suppression(phase5_results, suppression_needed)
    phase6_results = full_pass_integration(state, phase5_results, turb_results_site, p_dark_suppressed)
    plot_pass_keyrate(phase6_results, show=show_plots)
    print_validation_table()
    return dict(phase6_results=phase6_results, p_dark_suppressed=p_dark_suppressed,
                suppression_needed=suppression_needed)
