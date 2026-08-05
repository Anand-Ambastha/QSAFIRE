"""Part XI orchestrator: Final Validation Consolidation & Corrected Annual
Performance (Phase 8).

Reproduces notebook cells 213-224 in order. Consumes ``state`` (Part III),
``atm_state`` (Part IV), ``turb_state`` (Part V), ``phase5_state`` (Part VII),
``phase6_state`` (Part VIII), ``weather_state`` (Part IX), ``phase10_state``
(Part X). Returns a ``phase11_state`` dict.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.satellite.geometry import SITES, slant_from_elev
from src.atmosphere.transmittance import kasten_young_airmass, tau_aerosol_zenith, TAU_GAS_ZENITH
from src.atmosphere.turbulence import sigma_R2_slant, cn2_hv, LAMBDA_PHASE4_M
from src.atmosphere.seasonal import SEASONAL_AEROSOL, SEASONAL_TURBULENCE
from src.pat.losses import (
    ETA_RX_INTERNAL, ETA_DET_REAL, E_A_MISALIGNMENT, R_TX_HZ, POINTING_CASES,
    eta_diffraction, eta_pointing,
)
from src.optimization.joint import optimize_joint_ensemble

COLORS = {"Delhi (Alice)": "#1f77b4", "Mumbai (Bob)": "#d62728"}
SEASON_MAP_V2 = {"Winter": "Winter (DJF)", "Summer": "Summer/pre-monsoon (MAMJ)", "Monsoon": "Monsoon (JAS)"}


def recheck_seasonal_robustness(state, atm_state, phase5_results, turb_results_site, p_dark_suppressed,
                                 season_scale):
    """Cell 215: repeat Phase 7 Section 9.4's check with optimize_joint_ensemble."""
    results, idx = state['results'], state['idx']
    tau_R_810 = atm_state['tau_R_810']
    eta_point_val = eta_pointing(POINTING_CASES["Conservative default ([S2025])"]['sigma_j'])

    season_scale_v2 = {}
    for station in SITES:
        eta_best_ref = phase5_results[station]['eta_total'][idx].max()
        sigma_ref = float(turb_results_site[station]['sigma_R2'][idx].max())
        p_dark_ref = p_dark_suppressed[station]
        R_ref, _, _, _ = optimize_joint_ensemble(eta_best_ref, E_A_MISALIGNMENT, sigma_ref, p_dark_ref,
                                                   n_mc=15000, n_grid=9)

        L_m_best = results[station]['rng'][idx].max()*1000.0
        for season_short, season_full in SEASON_MAP_V2.items():
            aero = SEASONAL_AEROSOL[station][season_full]
            turb = SEASONAL_TURBULENCE[station][season_full]
            zenith_best = 90.0 - results[station]['el'][idx].max()
            airmass_best = kasten_young_airmass(zenith_best)
            tau_a = tau_aerosol_zenith(aero['aod0'], aero['alpha'], aero['lam0_nm'])
            tau_zenith_season = tau_R_810 + tau_a + TAU_GAS_ZENITH
            eta_atm_season = np.exp(-tau_zenith_season*airmass_best)
            eta_season = eta_diffraction(L_m_best) * eta_atm_season * eta_point_val * ETA_RX_INTERNAL * ETA_DET_REAL
            sigma_season = sigma_R2_slant(zenith_best, LAMBDA_PHASE4_M, lambda h, A0=turb['A0']: cn2_hv(h, A0=A0))

            R_season, _, _, _ = optimize_joint_ensemble(eta_season, E_A_MISALIGNMENT, sigma_season, p_dark_ref,
                                                          n_mc=15000, n_grid=9)
            ratio = max(R_season, 0.0) / R_ref if R_ref > 0 else 0.0
            season_scale_v2[(station, season_short)] = ratio

    print(f"{'Station':16s} {'Season':10s} {'OLD scale':>10s} {'NEW scale':>10s}  Status change")
    for station in SITES:
        for season in ["Winter", "Summer", "Monsoon"]:
            old_v = season_scale.get((station, season), 0.0)
            new_v = season_scale_v2.get((station, season), 0.0)
            change = "same (still zero)" if old_v == 0 and new_v == 0 else \
                     "NOW VIABLE" if old_v == 0 and new_v > 0 else \
                     "scaled" if old_v > 0 else "unexpected"
            print(f"{station:16s} {season:10s} {old_v:10.3f} {new_v:10.3f}  {change}")

    return season_scale_v2


def corrected_annual_key_volume(weather_state, phase10_results, season_scale_v2,
                                 save_path="phase11_annual_keyvolume_corrected.png", show=True):
    """Cells 218-219."""
    gated = weather_state['gated']
    pass_months, joint_avail, season_of_month = gated['pass_months'], gated['joint_avail'], gated['season_of_month']
    dt = None  # unused placeholder, kept for cell-order fidelity

    REFERENCE_BITS_PER_PASS_V2 = {}
    for station in SITES:
        R_clip = np.clip(phase10_results[station]['R'], 0, None)
        t = phase10_results[station]['t']
        dt_local = t[1] - t[0] if len(t) > 1 else 0.0
        REFERENCE_BITS_PER_PASS_V2[station] = np.sum(R_clip)*R_TX_HZ*dt_local

    annual_rows_v2 = []
    for station in SITES:
        total_bits_year = 0.0
        for season in ["Winter", "Summer", "Monsoon"]:
            mask = np.array([season_of_month[m] == season for m in pass_months])
            clear_passes = float(joint_avail[mask].sum()) if mask.sum() else 0.0
            scale = season_scale_v2.get((station, season), 0.0)
            bits_this_season = clear_passes * REFERENCE_BITS_PER_PASS_V2[station] * scale
            total_bits_year += bits_this_season
            annual_rows_v2.append(dict(station=station, season=season, expected_clear_passes=round(clear_passes, 1),
                                        season_scale=round(scale, 3), expected_key_bits=f"{bits_this_season:.2e}"))
        annual_rows_v2.append(dict(station=station, season="TOTAL/year", expected_clear_passes="-",
                                    season_scale="-", expected_key_bits=f"{total_bits_year:.2e}"))

    df_annual_v2 = pd.DataFrame(annual_rows_v2)
    print(df_annual_v2)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    seasons = ["Winter", "Summer", "Monsoon"]
    width = 0.35
    for k, station in enumerate(SITES):
        vals = []
        for season in seasons:
            mask = np.array([season_of_month[m] == season for m in pass_months])
            clear_passes = float(joint_avail[mask].sum()) if mask.sum() else 0.0
            scale = season_scale_v2.get((station, season), 0.0)
            vals.append(clear_passes * REFERENCE_BITS_PER_PASS_V2[station] * scale)
        ax.bar(np.arange(len(seasons)) + k*width, vals, width=width, label=station, color=COLORS[station])
    ax.set_xticks(np.arange(len(seasons)) + width/2)
    ax.set_xticklabels(seasons)
    ax.set_ylabel('Expected key bits / season / year')
    ax.set_title('Annual expected key volume by season -- corrected (joint signal+decoy optimization)')
    ax.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)

    return df_annual_v2, REFERENCE_BITS_PER_PASS_V2


def publication_keyrate_vs_elevation(state, atm_state, p_dark_suppressed,
                                      save_path="phase11_keyrate_vs_elevation.png", show=True):
    """Cell 221."""
    SITE_TAU_ZENITH = atm_state['SITE_TAU_ZENITH']
    eta_point_val = eta_pointing(POINTING_CASES["Conservative default ([S2025])"]['sigma_j'])

    el_sweep = np.linspace(20, 60, 25)
    R_vs_elev = {station: [] for station in SITES}
    for station in SITES:
        for el_test in el_sweep:
            zen = 90 - el_test
            air = kasten_young_airmass(zen)
            L_m_test = slant_from_elev(el_test, 500.0)*1000.0
            eta_test = (eta_diffraction(L_m_test) * np.exp(-SITE_TAU_ZENITH[station]*air)
                        * eta_point_val * ETA_RX_INTERNAL * ETA_DET_REAL)
            sigma_test = sigma_R2_slant(zen, LAMBDA_PHASE4_M, cn2_hv)
            R_test, _, _, _ = optimize_joint_ensemble(eta_test, E_A_MISALIGNMENT, sigma_test,
                                                        p_dark_suppressed[station], n_mc=10000, n_grid=8)
            R_vs_elev[station].append(max(R_test, 0.0)*R_TX_HZ/1000)

    fig, ax = plt.subplots(figsize=(8, 5))
    for station in SITES:
        ax.plot(el_sweep, R_vs_elev[station], color=COLORS[station], lw=2.2, marker='o', ms=3, label=station)
    ax.set_xlabel('Elevation angle (deg)')
    ax.set_ylabel('Key rate (kbit/s)')
    ax.set_title('Key rate vs. elevation angle (corrected, joint signal+decoy optimization)')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return R_vs_elev


def master_validation_table():
    """Cell 223."""
    master_validation = [
        dict(phase="Phase 2 (Part III)", check="Sun-sync inclination solver vs. Landsat-8/9, Sentinel-2",
             result="< 0.1 deg deviation", status="PASS"),
        dict(phase="Phase 2 (Part III)", check="Elevation<->slant-range formula vs. Liao et al. 2017 (Micius)",
             result="< 1.2% deviation", status="PASS"),
        dict(phase="Phase 2 (Part III)", check="Per-station pass duration vs. Micius reported ~273 s/day",
             result="282.5 s (Delhi) / 279.5 s (Mumbai)", status="PASS"),
        dict(phase="Phase 3 (Part IV)", check="Rayleigh zenith optical depth vs. AMS Glossary of Meteorology",
             result="0.1434 vs. 0.145 @ 500 nm (<2%)", status="PASS"),
        dict(phase="Phase 4 (Part V)", check="HV 5/7 canonical zenith Fried parameter",
             result="4.96 cm vs. 5.00 cm target (<1%)", status="PASS"),
        dict(phase="Phase 4 (Part V)", check="General slant Rytov-variance formula reduces to Part II horizontal formula",
             result="< 0.3% numerical-integration difference", status="PASS"),
        dict(phase="Phase 4 (Part V)", check="Delhi turbulence-classification correction (background dominance test)",
             result="84.8% of turbulence integral from below 1 km altitude", status="RESOLVED (see Section 5.8)"),
        dict(phase="Phase 5 (Part VII)", check="Diffraction loss formula vs. Liao et al. 2017 (~22 dB @ 1,200 km)",
             result="~19-31 dB over the pass's slant-range span; ~25 dB @ 1,200 km", status="PASS"),
        dict(phase="Phase 5 (Part VII)", check="Root-cause diagnostic: background light isolated as the binding constraint",
             result="Dark-count-only channel closes; night background alone reopens it", status="KEY FINDING"),
        dict(phase="Phase 6 (Part VIII)", check="Atmosphere+turbulence loss vs. Liao et al. 2017 (3-8 dB @ 1,200 km)",
             result="Delhi 3.03 dB, Mumbai 1.42 dB (shorter, better-case ranges)", status="PASS"),
        dict(phase="Phase 6 (Part VIII)", check="Per-pass key rate vs. real portable-ground-station literature",
             result="Matches 'order of 100 bps at 1,000 km' [Z2023] before Part X correction", status="PASS"),
        dict(phase="Phase 7 (Part IX)", check="Orbital pass-frequency computed directly from the Part III propagator",
             result="252/5,551 orbits/year within the ~15.5 deg capture range", status="COMPUTED"),
        dict(phase="Part X", check="Refactor correctness: joint rate formula at fixed decoy ratio vs. original",
             result="Exact match, diff < 1e-12", status="PASS"),
        dict(phase="Part X / XI", check="Joint signal+decoy optimization impact on per-pass key volume",
             result="Delhi +162%, Mumbai +108%", status="KEY FINDING"),
        dict(phase="Part XI", check="Seasonal robustness re-check with corrected optimizer",
             result="Mumbai now viable in all 3 seasons (was: winter only); Delhi non-viable in all 3", status="UPDATED FINDING"),
    ]
    df_master = pd.DataFrame(master_validation)
    print(df_master)
    return df_master


def run_full_part11_pipeline(state, atm_state, turb_state, phase5_state, phase6_state, weather_state,
                              phase10_state, show_plots=True):
    """Runs notebook cells 213-224 in order. Returns phase11_state dict."""
    season_scale_v2 = recheck_seasonal_robustness(
        state, atm_state, phase5_state['phase5_results'], turb_state['turb_results_site'],
        phase6_state['p_dark_suppressed'], weather_state['season_scale'])
    df_annual_v2, REFERENCE_BITS_PER_PASS_V2 = corrected_annual_key_volume(
        weather_state, phase10_state['phase10_results'], season_scale_v2, show=show_plots)
    R_vs_elev = publication_keyrate_vs_elevation(state, atm_state, phase6_state['p_dark_suppressed'], show=show_plots)
    df_master = master_validation_table()
    return dict(season_scale_v2=season_scale_v2, df_annual_v2=df_annual_v2,
                REFERENCE_BITS_PER_PASS_V2=REFERENCE_BITS_PER_PASS_V2, R_vs_elev=R_vs_elev, df_master=df_master)
