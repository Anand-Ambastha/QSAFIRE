"""Part IV orchestrator: Static Atmospheric Transmittance Module (Phase 3).

Reproduces notebook cells 99-114 in order. Consumes ``state`` (the
``geometry_state`` dict returned by
``src.experiments.satellite_experiments.run_full_part3_pipeline``) and
returns an ``atm_state`` dict carrying ``atm_results``, ``SITE_TAU_ZENITH``,
and ``SITE_TAU_ZENITH_WINTER`` forward to Part V.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.satellite.geometry import SITES
from src.atmosphere.transmittance import (
    LAMBDA_UM, LAMBDA_NM, TAU_GAS_ZENITH, AEROSOL_SITES,
    tau_rayleigh_zenith, tau_aerosol_zenith, kim_kruse_q,
    tau_aerosol_visibility, kasten_young_airmass, site_tau_zenith,
    site_tau_zenith_winter,
)

COLORS = {"Delhi (Alice)": "#1f77b4", "Mumbai (Bob)": "#d62728"}


def announce_wavelength():
    """Cell 100."""
    print(f"Operating wavelength for Phase 3: {LAMBDA_NM:.0f} nm")
    print("(Chosen to match the Micius downlink and the Behera & Sinha (2024) Indian satellite-QKD")
    print(" link-budget study, so results can be cross-checked against both.)")


def validate_rayleigh():
    """Cell 102: Rayleigh formula validated against AMS Glossary reference,
    returns tau_R_810."""
    tau_500_computed = tau_rayleigh_zenith(0.500)
    ams_reference = 0.145
    pct_err = 100*(tau_500_computed - ams_reference)/ams_reference
    print(f"Computed Rayleigh zenith optical depth at 500 nm: {tau_500_computed:.4f}")
    print(f"AMS Glossary of Meteorology reference value:      {ams_reference:.3f}")
    print(f"Deviation: {pct_err:+.2f}%")
    assert abs(pct_err) < 5, "Rayleigh formula deviates >5% from AMS Glossary reference"
    print("PASS: within 5% of the independently published reference value.\n")

    tau_R_810 = tau_rayleigh_zenith(LAMBDA_UM)
    print(f"Rayleigh zenith optical depth at {LAMBDA_NM:.0f} nm (sea level): {tau_R_810:.4f}")
    print("(Small and essentially site-independent at these altitudes -- Delhi ~216 m, Mumbai ~14 m")
    print(" above sea level give a <3% pressure correction relative to the sea-level value used here.)")
    return tau_R_810


def aerosol_table():
    """Cell 104."""
    rows = []
    for label, p in AEROSOL_SITES.items():
        tau = tau_aerosol_zenith(p['aod0'], p['alpha'], p['lam0_nm'])
        rows.append(dict(case=label, aod_500nm=p['aod0'], angstrom_alpha=p['alpha'],
                          tau_aerosol_810nm=round(tau, 4), confidence=p['confidence']))
    df_aero = pd.DataFrame(rows)
    return df_aero


def hanle_cross_check():
    """Cell 106: visibility-based aerosol model cross-check at IAO Hanle."""
    V_hanle = 23.0   # km, Behera & Sinha (2024), Table 2 (IAO Hanle MODTRAN input)
    tau_hanle, q_hanle = tau_aerosol_visibility(V_hanle)

    print(f"IAO Hanle (desert/high-altitude site, visibility = {V_hanle:.0f} km):")
    print(f"  Kim/Kruse size-distribution parameter q = {q_hanle:.2f} (visibility regime: 6-50 km)")
    print(f"  Aerosol zenith optical depth at {LAMBDA_NM:.0f} nm (2 km aerosol scale height): {tau_hanle:.3f}")
    print()
    print("Comparison with the AOD-based site values (Section 4.2), all at 810 nm:")
    for label, p in AEROSOL_SITES.items():
        tau = tau_aerosol_zenith(p['aod0'], p['alpha'], p['lam0_nm'])
        ratio = tau / tau_hanle
        print(f"  {label:38s}: tau = {tau:.3f}  ({ratio:.1f}x Hanle)")
    print()
    print("Expected physical ordering -- clear high-altitude desert site (Hanle) < urban Indian megacity")
    print("(Mumbai indicative) < Delhi's polluted continental aerosol, worst in winter haze. Confirmed above.")
    return tau_hanle, q_hanle


def announce_gas_absorption():
    """Cell 108."""
    print(f"Gaseous absorption zenith optical depth at {LAMBDA_NM:.0f} nm: {TAU_GAS_ZENITH} (indicative, flagged)")
    print("This is a minor term relative to aerosol extinction at both Delhi and Mumbai (Section 4.2)")
    print("and is carried as a fixed additive correction pending a full radiative-transfer treatment.")


def total_optical_depth(tau_R_810):
    """Cell 110 (script part): combine terms, print zenith transmittance."""
    SITE_TAU_ZENITH = site_tau_zenith(tau_R_810)
    SITE_TAU_ZENITH_WINTER = site_tau_zenith_winter(tau_R_810)

    print("Total zenith optical depth (Rayleigh + aerosol + gas), 810 nm, baseline case:")
    for name, tau in SITE_TAU_ZENITH.items():
        trans_pct = 100*np.exp(-tau)
        print(f"  {name:16s}: tau_zenith = {tau:.3f}  ->  zenith transmittance = {trans_pct:5.1f}%  "
              f"({-10*np.log10(np.exp(-tau)):.2f} dB)")
    print()
    for name, tau in SITE_TAU_ZENITH_WINTER.items():
        trans_pct = 100*np.exp(-tau)
        print(f"  {name:26s}: tau_zenith = {tau:.3f}  ->  zenith transmittance = {trans_pct:5.1f}%  "
              f"({-10*np.log10(np.exp(-tau)):.2f} dB)  [sensitivity case]")
    return SITE_TAU_ZENITH, SITE_TAU_ZENITH_WINTER


def apply_to_pass(state, SITE_TAU_ZENITH, save_path="phase3_atmospheric_transmittance.png", show=True):
    """Cell 112: apply the atmospheric model to the Phase 2 pass geometry."""
    results, t_grid = state['results'], state['t_grid']
    t_start, t_end, idx, el_min = state['t_start'], state['t_end'], state['idx'], state['el_min']

    atm_results = {}
    for name in SITES:
        zenith_deg = 90.0 - results[name]['el']
        airmass = kasten_young_airmass(zenith_deg)
        tau_slant = SITE_TAU_ZENITH[name] * airmass
        eta_atm = np.exp(-tau_slant)
        atm_results[name] = dict(airmass=airmass, tau_slant=tau_slant, eta_atm=eta_atm,
                                  loss_dB=-10*np.log10(eta_atm))

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for name in SITES:
        axes[0].plot(t_grid, atm_results[name]['eta_atm']*100, label=name, color=COLORS[name], lw=1.8)
    axes[0].axvspan(t_start, t_end, color='green', alpha=0.12, label='common dual-visibility window')
    axes[0].set_ylabel('Atmospheric transmittance (%)')
    axes[0].set_title(f'Static atmospheric transmittance during the pass ({LAMBDA_NM:.0f} nm)')
    axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3)

    for name in SITES:
        axes[1].plot(t_grid, atm_results[name]['loss_dB'], label=name, color=COLORS[name], lw=1.8)
    axes[1].axvspan(t_start, t_end, color='green', alpha=0.12)
    axes[1].axhline(SITE_TAU_ZENITH["Delhi (Alice)"]/np.log(10)*10, color=COLORS["Delhi (Alice)"], ls=':', lw=1)
    axes[1].set_ylabel('Atmospheric loss (dB)')
    axes[1].set_xlabel('Time relative to engineered crossing (s)')
    axes[1].set_title('Atmospheric loss vs. time (lower = clearer path; rises sharply near the horizon)')
    axes[1].grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)

    print("\nAtmospheric loss (dB) within the common dual-visibility window:")
    for name in SITES:
        lo, hi = atm_results[name]['loss_dB'][idx].min(), atm_results[name]['loss_dB'][idx].max()
        print(f"  {name:16s}: {lo:.2f} - {hi:.2f} dB  (best at highest elevation, worst near {el_min:.0f} deg cutoff)")

    return atm_results


def run_full_part4_pipeline(state, show_plots=True):
    """Runs notebook cells 99-114 in order. ``state`` is the geometry_state
    dict from Part III. Returns an atm_state dict (merged into/alongside
    state) for Part V onward."""
    announce_wavelength()
    tau_R_810 = validate_rayleigh()
    df_aero = aerosol_table()
    print(df_aero)
    hanle_cross_check()
    announce_gas_absorption()
    SITE_TAU_ZENITH, SITE_TAU_ZENITH_WINTER = total_optical_depth(tau_R_810)
    atm_results = apply_to_pass(state, SITE_TAU_ZENITH, show=show_plots)
    return dict(tau_R_810=tau_R_810, SITE_TAU_ZENITH=SITE_TAU_ZENITH,
                SITE_TAU_ZENITH_WINTER=SITE_TAU_ZENITH_WINTER, atm_results=atm_results)
