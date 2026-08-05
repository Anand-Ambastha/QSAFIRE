"""Part IX orchestrator: Weather-Gated Availability & Annual Performance (Phase 7).

Reproduces notebook cells 178-195 in order. Consumes ``state`` (Part III),
``atm_state`` (Part IV: SITE_TAU_ZENITH, tau_R_810), ``turb_state`` (Part V:
turb_results_site), ``phase5_state`` (Part VII: phase5_results),
``phase6_state`` (Part VIII: phase6_results, p_dark_suppressed). Returns a
``weather_state`` dict.
"""

import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.satellite.geometry import SITES, sub_satellite, slant_from_elev
from src.atmosphere.transmittance import kasten_young_airmass, tau_aerosol_zenith, TAU_GAS_ZENITH
from src.atmosphere.turbulence import sigma_R2_slant, cn2_hv, LAMBDA_PHASE4_M
from src.atmosphere.seasonal import SEASONAL_AEROSOL, SEASONAL_TURBULENCE, A0_BASELINE
from src.pat.losses import (
    LAMBDA_P5_M, OBSCURATION_RATIO, ETA_RX_INTERNAL, ETA_DET_REAL, E_A_MISALIGNMENT,
    R_TX_HZ, POINTING_CASES, eta_diffraction, eta_pointing,
)
from src.coupling.phase6 import optimize_ensemble_real
from src.weather.availability import RAINY_DAYS_PER_MONTH, FOG_DAYS_DELHI, DAYS_IN_MONTH, availability_fraction

COLORS = {"Delhi (Alice)": "#1f77b4", "Mumbai (Bob)": "#d62728"}

_offset_samples_deg = np.array([-9, -8, -7, -5, -3, -2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 7.5, 8, 9])
_duration_samples_s = np.array([0, 0, 46.0, 104.0, 132.5, 141.0, 146.0, 147.0, 147.5, 147.0, 146.5,
                                 144.5, 141.5, 138.0, 133.5, 121.5, 104.5, 79.5, 28.5, 14.0, 0, 0])


def duration_lookup(offset_deg):
    """Cell 182."""
    return np.clip(np.interp(offset_deg, _offset_samples_deg, _duration_samples_s,
                              left=0.0, right=0.0), 0.0, None)


def monthly_availability_table():
    """Cell 180 (script part)."""
    print("Monthly availability fraction (1 - outage proxy), by station:\n")
    print(f"{'Month':6s} {'Delhi':>8s} {'Mumbai':>8s}")
    for m in range(1, 13):
        print(f"{m:6d} {availability_fraction('Delhi (Alice)', m):8.3f} {availability_fraction('Mumbai (Bob)', m):8.3f}")


def orbital_pass_frequency(state):
    """Cell 182 (script part): full-year orbit-by-orbit longitude-offset
    distribution using the validated Part III propagator."""
    orb, a, inc, raan0, u0, gmst0 = state['orb'], state['a'], state['inc'], state['raan0'], state['u0'], state['gmst0']
    lat_m_unused = None
    # great-circle midpoint longitude, recomputed identically to Part III construct_pass()
    from src.satellite.geometry import gc_midpoint
    lat_m, lon_m = gc_midpoint(SITES["Delhi (Alice)"]['lat'], SITES["Delhi (Alice)"]['lon'],
                                SITES["Mumbai (Bob)"]['lat'],  SITES["Mumbai (Bob)"]['lon'])

    T_orbit = 2*np.pi/orb['u_dot']
    N_YEAR = int(365.25*86400/T_orbit)
    n_arr = np.arange(N_YEAR)
    t_arr = n_arr*T_orbit
    _, lon_c = sub_satellite(t_arr, a, inc, raan0, u0, orb['raan_dot'], orb['u_dot'], gmst0)
    dlon_year = (lon_c - lon_m + 180) % 360 - 180

    durations_year = duration_lookup(dlon_year)
    usable = durations_year > 0
    print(f"Orbital period (nodal): {T_orbit/60:.2f} min  ->  {N_YEAR} orbits/year")
    print(f"Orbits/year with a usable (duration>0) common-visibility pass: {usable.sum()}  "
          f"({100*usable.sum()/N_YEAR:.2f}% of all orbits)")
    print(f"Total annual common-visibility seconds (orbital geometry only, before weather gating): "
          f"{durations_year.sum():.0f} s  ({durations_year.sum()/usable.sum():.1f} s average per usable pass)")

    return dict(T_orbit=T_orbit, N_YEAR=N_YEAR, t_arr=t_arr, lon_m=lon_m,
                durations_year=durations_year, usable=usable)


def weather_gate_passes(orbit_freq, epoch=(2026, 7, 25, 12, 0, 0)):
    """Cell 184."""
    t_arr, usable, durations_year = orbit_freq['t_arr'], orbit_freq['usable'], orbit_freq['durations_year']
    epoch_dt = datetime.datetime(*epoch)
    pass_months = np.array([(epoch_dt + datetime.timedelta(seconds=float(t))).month for t in t_arr[usable]])
    pass_durations = durations_year[usable]

    avail_d = np.array([availability_fraction("Delhi (Alice)", m) for m in pass_months])
    avail_m = np.array([availability_fraction("Mumbai (Bob)", m) for m in pass_months])
    joint_avail = avail_d * avail_m  # independence assumption, Section 9.7

    expected_clear_seconds_per_pass = pass_durations * joint_avail
    total_expected_clear_seconds = expected_clear_seconds_per_pass.sum()
    effective_clear_passes = joint_avail.sum()

    print(f"Weather-gated annual result:")
    print(f"  Orbitally-usable passes/year: {usable.sum()}")
    print(f"  Expected fully-clear passes/year (weather-gated): {effective_clear_passes:.1f}")
    print(f"  Expected total clear common-visibility seconds/year: {total_expected_clear_seconds:.0f} s")
    print()
    print("By season (orbitally-usable passes, mean joint availability):")
    season_of_month = {12: "Winter", 1: "Winter", 2: "Winter", 3: "Summer", 4: "Summer", 5: "Summer",
                        6: "Monsoon", 7: "Monsoon", 8: "Monsoon", 9: "Monsoon", 10: "Post-monsoon", 11: "Post-monsoon"}
    for season in ["Winter", "Summer", "Monsoon", "Post-monsoon"]:
        mask = np.array([season_of_month[m] == season for m in pass_months])
        if mask.sum() == 0:
            continue
        print(f"  {season:14s}: {mask.sum():3d} orbitally-usable passes, "
              f"mean joint availability = {joint_avail[mask].mean():.3f}, "
              f"expected clear seconds = {expected_clear_seconds_per_pass[mask].sum():.0f} s")

    return dict(pass_months=pass_months, pass_durations=pass_durations, joint_avail=joint_avail,
                season_of_month=season_of_month)


def per_pass_key_volume_by_season(state, phase5_results, turb_results_site, phase6_results,
                                   p_dark_suppressed, tau_R_810, SITE_TAU_ZENITH):
    """Cells 186 (script part): per-pass key volume by season, rescaled from
    the Phase 6 reference case."""
    results, idx = state['results'], state['idx']
    dt = state['t_grid'][1] - state['t_grid'][0]
    eta_point_val = eta_pointing(POINTING_CASES["Conservative default ([S2025])"]['sigma_j'])

    REFERENCE_BITS_PER_PASS = {station: np.sum(np.clip(phase6_results[station]['R'], 0, None))*R_TX_HZ*dt
                                for station in SITES}
    print("Phase 6 reference per-pass bits (baseline case):", REFERENCE_BITS_PER_PASS)

    SEASON_MAP_PHASE6 = {"Winter": "Winter (DJF)", "Summer": "Summer/pre-monsoon (MAMJ)", "Monsoon": "Monsoon (JAS)"}

    season_scale = {}
    for station in SITES:
        eta_best_ref = phase5_results[station]['eta_total'][idx].max()
        sigma_ref = float(turb_results_site[station]['sigma_R2'][idx].max())
        p_dark_ref = p_dark_suppressed[station]
        R_ref, _, _ = optimize_ensemble_real(eta_best_ref, E_A_MISALIGNMENT, sigma_ref, p_dark_ref, n_mc=20000, n_grid=12)

        L_m_best = results[station]['rng'][idx].max()*1000.0
        for season_short, season_full in SEASON_MAP_PHASE6.items():
            aero = SEASONAL_AEROSOL[station][season_full]
            turb = SEASONAL_TURBULENCE[station][season_full]
            zenith_best = 90.0 - results[station]['el'][idx].max()
            airmass_best = kasten_young_airmass(zenith_best)
            tau_a = tau_aerosol_zenith(aero['aod0'], aero['alpha'], aero['lam0_nm'])
            tau_zenith_season = tau_R_810 + tau_a + TAU_GAS_ZENITH
            eta_atm_season = np.exp(-tau_zenith_season*airmass_best)
            eta_season = eta_diffraction(L_m_best) * eta_atm_season * eta_point_val * ETA_RX_INTERNAL * ETA_DET_REAL
            sigma_season = sigma_R2_slant(zenith_best, LAMBDA_PHASE4_M, lambda h, A0=turb['A0']: cn2_hv(h, A0=A0))

            R_season, _, _ = optimize_ensemble_real(eta_season, E_A_MISALIGNMENT, sigma_season, p_dark_ref,
                                                      n_mc=20000, n_grid=12)
            ratio = max(R_season, 0.0) / R_ref if R_ref > 0 else 0.0
            season_scale[(station, season_short)] = ratio
            print(f"{station:16s} {season_short:8s}: R_season/R_reference = {ratio:.3f}")

    return REFERENCE_BITS_PER_PASS, season_scale


def annual_expected_key_volume(gated, REFERENCE_BITS_PER_PASS, season_scale,
                                save_path="phase7_annual_keyvolume.png", show=True):
    """Cells 189-190."""
    pass_months, joint_avail, season_of_month = gated['pass_months'], gated['joint_avail'], gated['season_of_month']

    annual_rows = []
    for station in SITES:
        total_bits_year = 0.0
        for season in ["Winter", "Summer", "Monsoon"]:
            mask = np.array([season_of_month[m] == season for m in pass_months])
            n_orbital_passes = int(mask.sum())
            clear_passes = float(joint_avail[mask].sum()) if mask.sum() else 0.0
            scale = season_scale.get((station, season), 0.0)
            bits_this_season = clear_passes * REFERENCE_BITS_PER_PASS[station] * scale
            total_bits_year += bits_this_season
            annual_rows.append(dict(station=station, season=season, orbital_passes=n_orbital_passes,
                                     expected_clear_passes=round(clear_passes, 1),
                                     season_scale=round(scale, 3),
                                     expected_key_bits=f"{bits_this_season:.2e}"))
        annual_rows.append(dict(station=station, season="TOTAL/year", orbital_passes="-", expected_clear_passes="-",
                                 season_scale="-", expected_key_bits=f"{total_bits_year:.2e}"))

    df_annual = pd.DataFrame(annual_rows)
    print(df_annual)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for station in SITES:
        seasons = ["Winter", "Summer", "Monsoon"]
        vals = []
        for season in seasons:
            mask = np.array([season_of_month[m] == season for m in pass_months])
            clear_passes = float(joint_avail[mask].sum()) if mask.sum() else 0.0
            scale = season_scale.get((station, season), 0.0)
            vals.append(clear_passes * REFERENCE_BITS_PER_PASS[station] * scale)
        ax.bar(np.arange(len(seasons)) + (0 if "Delhi" in station else 0.35), vals, width=0.35,
               label=station, color=COLORS[station])
    ax.set_xticks(np.arange(len(seasons)) + 0.175)
    ax.set_xticklabels(seasons)
    ax.set_ylabel('Expected key bits / season / year')
    ax.set_title('Annual expected key volume by season (weather-gated)')
    ax.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)

    return df_annual


def sensitivity_analysis(state, phase5_results, turb_results_site, p_dark_suppressed, atm_results,
                          SITE_TAU_ZENITH, save_path="phase7_sensitivity_tornado.png", show=True):
    """Cells 192-193: one-station, one-point sensitivity sweep across three
    design levers (min elevation, receiver aperture, ground turbulence A0)."""
    results, idx = state['results'], state['idx']
    eta_point_val = eta_pointing(POINTING_CASES["Conservative default ([S2025])"]['sigma_j'])
    station = "Mumbai (Bob)"  # representative station for the sweep
    eta_best_ref = phase5_results[station]['eta_total'][idx].max()
    sigma_ref = float(turb_results_site[station]['sigma_R2'][idx].max())
    p_dark_ref = p_dark_suppressed[station]
    R_ref, _, _ = optimize_ensemble_real(eta_best_ref, E_A_MISALIGNMENT, sigma_ref, p_dark_ref, n_mc=20000, n_grid=12)

    L_m_best = results[station]['rng'][idx].max()*1000.0
    sensitivity_rows = []

    # --- Lever 1: minimum elevation angle ---
    for el_test in [15, 20, 25, 30]:
        zen = 90 - el_test
        air = kasten_young_airmass(zen)
        L_m_test = slant_from_elev(el_test, 500.0)*1000.0
        eta_test = eta_diffraction(L_m_test) * np.exp(-SITE_TAU_ZENITH[station]*air) * eta_point_val * ETA_RX_INTERNAL * ETA_DET_REAL
        R_test, _, _ = optimize_ensemble_real(eta_test, E_A_MISALIGNMENT, sigma_ref, p_dark_ref, n_mc=20000, n_grid=12)
        sensitivity_rows.append(dict(lever="min elevation angle", value=f"{el_test} deg",
                                      R=R_test, pct_change=100*(R_test-R_ref)/R_ref if R_ref else np.nan))

    # --- Lever 2: receiver aperture D_RX ---
    for D_test in [0.7, 1.0, 1.3, 1.6]:
        D_int_test = D_test*OBSCURATION_RATIO
        G_rx_test = (np.pi**2/LAMBDA_P5_M**2)*(D_test**2 - D_int_test**2)
        eta_diff_test = eta_diffraction(L_m_best, G_rx=G_rx_test)
        eta_test = eta_diff_test * atm_results[station]['eta_atm'][idx].max() * eta_point_val * ETA_RX_INTERNAL * ETA_DET_REAL
        R_test, _, _ = optimize_ensemble_real(eta_test, E_A_MISALIGNMENT, sigma_ref, p_dark_ref, n_mc=20000, n_grid=12)
        sensitivity_rows.append(dict(lever="receiver aperture D_Rx", value=f"{D_test} m",
                                      R=R_test, pct_change=100*(R_test-R_ref)/R_ref if R_ref else np.nan))

    # --- Lever 3: ground-level turbulence A0 ---
    for A0_test in [0.5*A0_BASELINE, A0_BASELINE, 2*A0_BASELINE, 7*A0_BASELINE]:
        sigma_test = sigma_R2_slant(90-results[station]['el'][idx].max(), LAMBDA_PHASE4_M,
                                     lambda h, A0=A0_test: cn2_hv(h, A0=A0))
        R_test, _, _ = optimize_ensemble_real(eta_best_ref, E_A_MISALIGNMENT, sigma_test, p_dark_ref, n_mc=20000, n_grid=12)
        sensitivity_rows.append(dict(lever="ground turbulence A0", value=f"{A0_test:.1e}",
                                      R=R_test, pct_change=100*(R_test-R_ref)/R_ref if R_ref else np.nan))

    df_sens = pd.DataFrame(sensitivity_rows)
    print(df_sens)

    fig, ax = plt.subplots(figsize=(8, 5))
    lever_ranges = df_sens.groupby('lever')['pct_change'].agg(lambda x: x.max()-x.min()).sort_values()
    ax.barh(lever_ranges.index, lever_ranges.values, color=['#2a6f97', '#d62728', '#2ca02c'])
    ax.set_xlabel('Range of key-rate change across the tested values (%)')
    ax.set_title(f'Sensitivity comparison at {station} (best-elevation point)')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)

    print("Most sensitive lever:", lever_ranges.index[-1], f"({lever_ranges.values[-1]:.0f}% range)")
    print("Least sensitive lever:", lever_ranges.index[0], f"({lever_ranges.values[0]:.0f}% range)")

    return df_sens


def run_full_part9_pipeline(state, atm_state, turb_state, phase5_state, phase6_state, show_plots=True):
    """Runs notebook cells 178-195 in order. Returns weather_state dict."""
    monthly_availability_table()
    orbit_freq = orbital_pass_frequency(state)
    gated = weather_gate_passes(orbit_freq)
    REFERENCE_BITS_PER_PASS, season_scale = per_pass_key_volume_by_season(
        state, phase5_state['phase5_results'], turb_state['turb_results_site'], phase6_state['phase6_results'],
        phase6_state['p_dark_suppressed'], atm_state['tau_R_810'], atm_state['SITE_TAU_ZENITH'])
    df_annual = annual_expected_key_volume(gated, REFERENCE_BITS_PER_PASS, season_scale, show=show_plots)
    df_sens = sensitivity_analysis(state, phase5_state['phase5_results'], turb_state['turb_results_site'],
                                    phase6_state['p_dark_suppressed'], atm_state['atm_results'],
                                    atm_state['SITE_TAU_ZENITH'], show=show_plots)
    return dict(orbit_freq=orbit_freq, gated=gated, REFERENCE_BITS_PER_PASS=REFERENCE_BITS_PER_PASS,
                season_scale=season_scale, df_annual=df_annual, df_sens=df_sens)
