"""Part V orchestrator: Slant-Path Turbulence & Scintillation Module (Phase 4).

Reproduces notebook cells 115-130 in order. Consumes ``state`` (geometry_state
from Part III) and returns a ``turb_state`` dict carrying ``turb_results``,
``turb_results_site``, ``gg_results``, ``gg_results_site`` forward to Part VI
onward.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.satellite.geometry import SITES
from src.atmosphere.turbulence import (
    LAMBDA_PHASE4_M, K_810, H_TOP_M, SPHERICAL_UPLINK_FACTOR, D_PROVISIONAL_M,
    TURBULENCE_SITE_CASES, cn2_hv, fried_r0, sigma_R2_slant_planewave,
    sigma_R2_slant, aoa_jitter_urad, strehl_uncompensated,
)
from src.fso.config import WAVELENGTH_M
from src.fso.channel_models import rytov_variance, gamma_gamma_params, classify_turbulence_regime

COLORS = {"Delhi (Alice)": "#1f77b4", "Mumbai (Bob)": "#d62728"}


def plot_hv_profile(save_path="phase4_hv_profile.png", show=True):
    """Cell 116 (script part)."""
    h_plot = np.linspace(0, H_TOP_M, 2000)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(cn2_hv(h_plot), h_plot/1000, color='#2a6f97')
    ax.set_xscale('log')
    ax.set_xlabel(r'$C_n^2(h)$  (m$^{-2/3}$)')
    ax.set_ylabel('Altitude (km)')
    ax.set_title('Hufnagel-Valley 5/7 turbulence profile (A0=1.7e-14, v=21 m/s)')
    ax.grid(alpha=0.3, which='both')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


def validate_fried_r0():
    """Cell 118: HV 5/7 canonical zenith Fried parameter validation."""
    r0_500_zenith, J_cn2 = fried_r0(0.0, 500e-9)
    print(f"Integrated turbulence strength  J = int Cn2 dh (0-30 km) = {J_cn2:.4e} m^(1/3)")
    print(f"Zenith Fried parameter r0 @ 500 nm: {r0_500_zenith*100:.2f} cm   (HV 5/7 canonical target: 5.00 cm)")
    pct_err = 100*(r0_500_zenith*100 - 5.0)/5.0
    print(f"Deviation: {pct_err:+.2f}%")
    assert abs(pct_err) < 5, "HV profile r0 deviates >5% from the canonical HV 5/7 target"
    print("PASS: within 5% of the canonical HV 5/7 reference value.")
    return r0_500_zenith, J_cn2


def validate_slant_rytov():
    """Cell 120 (script part): two-stage validation of the plane-wave
    machinery and the corrected spherical-wave (uplink) output."""
    _Cn2_const_test = 1e-14
    _L_test_km = 6.0
    _sigma_new_plane = sigma_R2_slant_planewave(0.0, WAVELENGTH_M, lambda h: _Cn2_const_test, h0=0.0, h_top=_L_test_km*1000.0)
    _sigma_old_plane = rytov_variance(_L_test_km, _Cn2_const_test)
    _pct_plane = 100*(_sigma_new_plane - _sigma_old_plane)/_sigma_old_plane
    print(f"[Stage 1: plane-wave path-integral machinery, unchanged]")
    print(f"General path-integral formula (constant Cn2, zenith):  sigma_R^2 = {_sigma_new_plane:.5f}")
    print(f"Existing Part II horizontal formula:                    sigma_R^2 = {_sigma_old_plane:.5f}")
    print(f"Relative difference: {_pct_plane:+.3f}%  (numerical-integration artifact only)")
    assert abs(_pct_plane) < 1.0, "Plane-wave path-integral machinery no longer reduces to the existing horizontal formula"
    print("PASS: plane-wave path-integral machinery unchanged and still correct.\n")

    _sigma_new_spherical = sigma_R2_slant(0.0, WAVELENGTH_M, lambda h: _Cn2_const_test, h0=0.0, h_top=_L_test_km*1000.0)
    _beta0_target = 0.5 * _Cn2_const_test * (2*np.pi/WAVELENGTH_M)**(7/6) * (_L_test_km*1000.0)**(11/6)
    _pct_spherical = 100*(_sigma_new_spherical - _beta0_target)/_beta0_target
    print(f"[Stage 2: corrected spherical-wave / uplink output]")
    print(f"Corrected sigma_R^2 (spherical/uplink):                 sigma_R^2 = {_sigma_new_spherical:.5f}")
    print(f"Textbook spherical-wave horizontal target (0.5 coeff):  sigma_R^2 = {_beta0_target:.5f}")
    print(f"Relative difference: {_pct_spherical:+.3f}%  (numerical-integration artifact only)")
    assert abs(_pct_spherical) < 1.0, "Corrected spherical-wave slant formula does not reduce to the textbook spherical-wave horizontal value"
    print("PASS: corrected general formula reduces to the textbook spherical-wave (uplink) horizontal value.")


def turbulence_site_cases_table():
    """Cell 122."""
    for label, p in TURBULENCE_SITE_CASES.items():
        r0z, _ = fried_r0(0.0, LAMBDA_PHASE4_M, A0=p['A0'], v=p['v'])
        print(f"{label:38s}: A0={p['A0']:.2e}  ->  zenith r0 @ 810 nm = {r0z*100:5.2f} cm   [{p['confidence']}]")


def apply_to_pass(state, save_path="phase4_turbulence_vs_time.png", show=True):
    """Cells 124: apply the turbulence model to the Phase 2 pass geometry,
    shared baseline vs. Delhi site case."""
    results, t_grid = state['results'], state['t_grid']
    t_start, t_end, idx = state['t_start'], state['t_end'], state['idx']

    SITE_A0 = {
        "Delhi (Alice)": TURBULENCE_SITE_CASES["Delhi urban heat-island (indicative)"]['A0'],
        "Mumbai (Bob)":  TURBULENCE_SITE_CASES["shared baseline (HV 5/7)"]['A0'],
    }

    turb_results = {}       # shared baseline for BOTH sites (literature default, no site data)
    turb_results_site = {}  # site-specific: Delhi enhanced, Mumbai baseline (Section 5.3 sensitivity case)
    for name in SITES:
        zenith_deg = 90.0 - results[name]['el']
        sigma_R2_t = np.array([sigma_R2_slant(z, LAMBDA_PHASE4_M, cn2_hv) for z in zenith_deg])
        r0_t = np.array([fried_r0(z, LAMBDA_PHASE4_M)[0] for z in zenith_deg])
        turb_results[name] = dict(sigma_R2=sigma_R2_t, r0=r0_t)

        A0_site = SITE_A0[name]
        sigma_R2_site = np.array([sigma_R2_slant(z, LAMBDA_PHASE4_M, lambda h, A0=A0_site: cn2_hv(h, A0=A0))
                                   for z in zenith_deg])
        r0_site = np.array([fried_r0(z, LAMBDA_PHASE4_M, A0=A0_site)[0] for z in zenith_deg])
        turb_results_site[name] = dict(sigma_R2=sigma_R2_site, r0=r0_site)

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for name in SITES:
        axes[0].plot(t_grid, turb_results[name]['sigma_R2'], color=COLORS[name], lw=1.2, ls='--',
                     label=f'{name} (shared baseline)')
        axes[0].plot(t_grid, turb_results_site[name]['sigma_R2'], color=COLORS[name], lw=2.0,
                     label=f'{name} (site case, A0={SITE_A0[name]:.1e})')
    axes[0].axhline(1.0, color='gray', ls=':', lw=1)
    axes[0].text(t_grid[0], 1.05, 'weak / moderate boundary', fontsize=8, color='gray')
    axes[0].axvspan(t_start, t_end, color='green', alpha=0.12, label='common dual-visibility window')
    axes[0].set_ylabel(r'Rytov variance $\sigma_R^2$')
    axes[0].set_title('Slant-path turbulence strength during the pass (810 nm): baseline vs. site case')
    axes[0].legend(fontsize=8, ncol=2); axes[0].grid(alpha=0.3)

    for name in SITES:
        axes[1].plot(t_grid, turb_results_site[name]['r0']*100, color=COLORS[name], lw=2.0, label=name)
    axes[1].axvspan(t_start, t_end, color='green', alpha=0.12)
    axes[1].set_ylabel(r'Fried parameter $r_0$ (cm)')
    axes[1].set_xlabel('Time relative to engineered crossing (s)')
    axes[1].set_title('Slant-path Fried parameter vs. time (site case)')
    axes[1].grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)

    print("Rytov variance in the common dual-visibility window -- shared baseline vs. site case:")
    for name in SITES:
        s2b = turb_results[name]['sigma_R2'][idx]
        s2s = turb_results_site[name]['sigma_R2'][idx]
        print(f"  {name:16s} baseline: [{s2b.min():.3f}, {s2b.max():.3f}]   "
              f"site case: [{s2s.min():.3f}, {s2s.max():.3f}]")
    print()
    print("Delhi's site case crosses sigma_R^2 = 1 (weak/moderate boundary) near the low-elevation edge")
    print("of the window -- i.e. a plausible (not yet measured) urban enhancement is already enough to")
    print("push Delhi into the moderate-turbulence regime for part of the pass, even though the shared,")
    print("unenhanced literature baseline alone would call the whole pass 'weak' for both cities.")

    return turb_results, turb_results_site, SITE_A0


def feed_gamma_gamma(state, turb_results, turb_results_site):
    """Cell 126: feed sigma_R^2(t) into the existing (Part II) Gamma-Gamma
    machinery, shared baseline and Delhi site case both retained."""
    idx = state['idx']
    gg_results = {}       # shared baseline
    gg_results_site = {}  # site case (Delhi enhanced, Mumbai baseline)
    for name in SITES:
        for store, source in [(gg_results, turb_results), (gg_results_site, turb_results_site)]:
            s2 = source[name]['sigma_R2']
            ag, bg = gamma_gamma_params(s2)
            regime = classify_turbulence_regime(s2)
            store[name] = dict(alpha_g=ag, beta_g=bg, regime=np.asarray(regime))

    for name in SITES:
        print(f"{name}:")
        for label, store in [("shared baseline", gg_results), ("site case", gg_results_site)]:
            regimes_in_window, counts = np.unique(store[name]['regime'][idx], return_counts=True)
            regime_summary = ", ".join(f"{r}: {100*c/len(idx):.0f}%" for r, c in zip(regimes_in_window, counts))
            ag_mean = store[name]['alpha_g'][idx].mean()
            bg_mean = store[name]['beta_g'][idx].mean()
            print(f"  [{label:15s}] regime -> {regime_summary:35s}  "
                  f"mean alpha_g={ag_mean:.2f}, beta_g={bg_mean:.2f}")

    return gg_results, gg_results_site


def provisional_pat(state, turb_results_site):
    """Cell 128 (script part): provisional AOA jitter / Strehl ratio."""
    results, idx = state['results'], state['idx']
    for name in SITES:
        zenith_deg = 90.0 - results[name]['el'][idx]
        aoa = np.array([aoa_jitter_urad(z) for z in zenith_deg])
        strehl = strehl_uncompensated(turb_results_site[name]['r0'][idx])
        print(f"{name:16s} (common-visibility window, D={D_PROVISIONAL_M:.1f} m provisional, site-case r0):")
        print(f"{'':16s}   AOA jitter        = [{aoa.min():.2f}, {aoa.max():.2f}] urad rms")
        print(f"{'':16s}   Uncompensated S   = [{strehl.min():.2e}, {strehl.max():.2e}]  "
              f"(D/r0 >> 1 -> strong AO or aperture-averaging correction required; carried to Phase 5)")


def run_full_part5_pipeline(state, show_plots=True):
    """Runs notebook cells 115-130 in order. Returns turb_state dict for
    Part VI onward."""
    plot_hv_profile(show=show_plots)
    validate_fried_r0()
    validate_slant_rytov()
    turbulence_site_cases_table()
    turb_results, turb_results_site, SITE_A0 = apply_to_pass(state, show=show_plots)
    gg_results, gg_results_site = feed_gamma_gamma(state, turb_results, turb_results_site)
    provisional_pat(state, turb_results_site)
    return dict(turb_results=turb_results, turb_results_site=turb_results_site,
                gg_results=gg_results, gg_results_site=gg_results_site, SITE_A0=SITE_A0)
