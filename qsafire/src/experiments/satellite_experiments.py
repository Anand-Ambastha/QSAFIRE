"""Part III orchestrator: Ground-to-satellite geometry engine (Phase 2).

Reproduces notebook cells 79-98 in order. Returns a ``geometry_state`` dict
that later parts (IV onward) consume, exactly mirroring the notebook's use of
module-level globals (``results``, ``site_ecef``, ``orb``, ``a``, ``inc``,
``raan0``, ``u0``, ``gmst0``, ``t_grid``, ``common``, ``idx``, ``t_start``,
``t_end``) across narrative sections.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.satellite.geometry import (
    MU_EARTH, RE_EQ, F_WGS84, E2_WGS84, J2, OMEGA_E, YEAR_TROPICAL_DAYS,
    SITES, geodetic_to_ecef, enu_rotation, build_site_ecef, sun_sync_orbit,
    julian_date, gmst_rad, eci_position, eci_to_ecef, look_angles,
    refraction_bump_deg, slant_from_elev, gc_midpoint, sub_satellite,
    haversine_km,
)

EL_MIN = 20.0   # baseline minimum operational elevation (Phase 1, Section 3.3)


def print_sites_table():
    """Cell 80 (script part)."""
    print(pd.DataFrame(SITES).T)
    return pd.DataFrame(SITES).T


def geodesy_ecef_check():
    """Cell 82 (script part)."""
    site_ecef = build_site_ecef()
    for name, r in site_ecef.items():
        print(f"{name:16s} ECEF = ({r[0]:9.2f}, {r[1]:9.2f}, {r[2]:9.2f}) km   |r| = {np.linalg.norm(r):.2f} km")
    print()
    print(f"Sanity check: |r| should equal ~Earth radius (6378 eq / 6357 polar). Both stations are near-")
    print(f"equatorial-ish latitudes so |r| should sit close to {RE_EQ:.1f} km.")
    return site_ecef


def orbit_design_table():
    """Cell 84 (script part): sun-sync solver, candidate + real-mission checks."""
    rows = []
    for label, h in [("candidate: 450 km", 450), ("candidate: 500 km (Micius altitude)", 500),
                     ("candidate: 550 km", 550), ("candidate: 600 km", 600),
                     ("check: Landsat-8/9 (705 km)", 705), ("check: Sentinel-2 (786 km)", 786)]:
        o = sun_sync_orbit(h)
        rows.append(dict(case=label, alt_km=h, inclination_deg=np.degrees(o['i']),
                          period_min=2*np.pi/o['n0']/60, raan_dot_deg_per_day=np.degrees(o['raan_dot'])*86400))
    df_orbit = pd.DataFrame(rows)
    return df_orbit


def validate_slant_range():
    """Cell 89: validate elevation<->slant-range formula against Liao et al. (2017)."""
    val_rows = []
    for elev, h, ref in [(85.7, 500.0, 507.0), (25.0, 500.0, 1034.7)]:
        d = slant_from_elev(elev, h)
        val_rows.append(dict(elevation_deg=elev, altitude_km=h, computed_km=round(d, 2),
                              liao_2017_reported_km=ref, pct_diff=round(100*(d-ref)/ref, 2)))
    df_val1 = pd.DataFrame(val_rows)
    print(df_val1.to_string(index=False))
    assert (df_val1['pct_diff'].abs() < 2.0).all(), "Elevation<->range formula deviates >2% from Micius benchmark"
    print()
    print("PASS: within 2% of Liao et al. (2017) published slant-range values.")
    return df_val1


def construct_pass(h_alt=500.0, epoch=(2026, 7, 25, 12, 0, 0), el_min=EL_MIN):
    """Cell 91 + 92: engineer a Delhi-Mumbai dual-downlink crossing at the
    great-circle midpoint, then propagate elevation/azimuth/slant-range for
    both stations across the pass. Returns the full ``geometry_state`` dict."""
    site_ecef = build_site_ecef()
    orb = sun_sync_orbit(h_alt)
    a, inc = orb['a'], orb['i']

    lat_m, lon_m = gc_midpoint(SITES["Delhi (Alice)"]['lat'], SITES["Delhi (Alice)"]['lon'],
                                SITES["Mumbai (Bob)"]['lat'],  SITES["Mumbai (Bob)"]['lon'])

    d_gc = haversine_km(SITES["Delhi (Alice)"]['lat'], SITES["Delhi (Alice)"]['lon'],
                         SITES["Mumbai (Bob)"]['lat'],  SITES["Mumbai (Bob)"]['lon'])
    print(f"Delhi-Mumbai great-circle baseline: {d_gc:.1f} km  (Phase 1 document, Section 3.3, quotes 1,148 km)")
    print(f"Great-circle midpoint: lat={lat_m:.3f}, lon={lon_m:.3f}")

    jd0 = julian_date(*epoch)
    gmst0 = gmst_rad(jd0)

    lat_t, lon_t = np.radians(lat_m), np.radians(lon_m)
    u0 = np.arcsin(np.sin(lat_t)/np.sin(inc))                                  # ascending branch
    lon_arg = np.arctan2(np.cos(inc)*np.sin(u0), np.cos(u0))
    raan0 = lon_t + gmst0 - lon_arg

    lat_chk, lon_chk = sub_satellite(0.0, a, inc, raan0, u0, orb['raan_dot'], orb['u_dot'], gmst0)
    print(f"Sub-satellite point at t=0 (should equal the midpoint above): lat={lat_chk:.3f}, lon={lon_chk:.3f}")

    # Propagate elevation / azimuth / slant range for both stations across the pass
    t_grid = np.linspace(-360, 360, 1441)  # seconds around the engineered crossing, 0.5 s steps

    results = {name: {'el': [], 'az': [], 'rng': []} for name in SITES}
    for t in t_grid:
        raan = raan0 + orb['raan_dot']*t
        u = u0 + orb['u_dot']*t
        r_eci = eci_position(a, inc, raan, u)
        gmst = gmst0 + OMEGA_E*t
        r_ecef = eci_to_ecef(r_eci, gmst)
        for name, s in SITES.items():
            el, az, rng = look_angles(r_ecef, site_ecef[name], s['lat'], s['lon'])
            results[name]['el'].append(el); results[name]['az'].append(az); results[name]['rng'].append(rng)
    results = {name: {k: np.array(v) for k, v in d.items()} for name, d in results.items()}

    common = np.ones_like(t_grid, dtype=bool)
    for name in SITES:
        common &= results[name]['el'] >= el_min

    dt = t_grid[1]-t_grid[0]
    print(f"Common dual-visibility window (elev >= {el_min:.0f} deg at BOTH stations): "
          f"{common.sum()} samples = {common.sum()*dt:.1f} s")
    idx = np.where(common)[0]
    t_start, t_end = t_grid[idx[0]], t_grid[idx[-1]]
    print(f"Window spans t = {t_start:+.1f} s to {t_end:+.1f} s relative to the engineered crossing.")
    print()

    summary_rows = []
    for name in SITES:
        own_mask = results[name]['el'] >= el_min
        elev_lo, elev_hi = results[name]['el'][idx].min(), results[name]['el'][idx].max()
        rng_lo, rng_hi = results[name]['rng'][idx].min(), results[name]['rng'][idx].max()
        summary_rows.append(dict(
            station=name,
            own_pass_duration_s=round(own_mask.sum()*dt, 1),
            common_window_duration_s=round(common.sum()*dt, 1),
            max_elevation_deg=round(results[name]['el'].max(), 2),
            elev_in_common_window_deg=f"[{elev_lo:.1f}, {elev_hi:.1f}]",
            slant_range_in_common_window_km=f"[{rng_lo:.1f}, {rng_hi:.1f}]",
        ))
    df_pass = pd.DataFrame(summary_rows)

    return dict(
        site_ecef=site_ecef, orb=orb, a=a, inc=inc, raan0=raan0, u0=u0, gmst0=gmst0,
        t_grid=t_grid, results=results, common=common, idx=idx, t_start=t_start, t_end=t_end,
        df_pass=df_pass, el_min=el_min, h_alt=h_alt,
    )


def plot_pass_geometry(state, save_path="phase2_pass_geometry.png", show=True):
    """Cell 94."""
    t_grid, results = state['t_grid'], state['results']
    t_start, t_end, el_min = state['t_start'], state['t_end'], state['el_min']
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    colors = {"Delhi (Alice)": "#1f77b4", "Mumbai (Bob)": "#d62728"}
    for name in SITES:
        axes[0].plot(t_grid, results[name]['el'], label=name, color=colors[name], lw=1.8)
    axes[0].axhline(el_min, color='gray', ls='--', lw=1, label=f'min. elevation ({el_min:.0f} deg)')
    axes[0].axvspan(t_start, t_end, color='green', alpha=0.12, label='common dual-visibility window')
    axes[0].set_ylabel('Elevation (deg)')
    axes[0].set_title('Delhi-Mumbai dual-downlink pass: elevation angle vs. time')
    axes[0].legend(loc='upper right', fontsize=9)
    axes[0].grid(alpha=0.3)

    for name in SITES:
        axes[1].plot(t_grid, results[name]['rng'], label=name, color=colors[name], lw=1.8)
    axes[1].axvspan(t_start, t_end, color='green', alpha=0.12)
    axes[1].set_ylabel('Slant range (km)')
    axes[1].set_xlabel('Time relative to engineered crossing (s)')
    axes[1].set_title('Slant range vs. time')
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


def day_sweep(state, save_path="phase2_day_sweep.png", show=True):
    """Cell 96: 24-hour elevation history for the same fixed orbital plane."""
    a, inc, raan0, u0, gmst0 = state['a'], state['inc'], state['raan0'], state['u0'], state['gmst0']
    orb, site_ecef, el_min = state['orb'], state['site_ecef'], state['el_min']

    t_day = np.linspace(-12*3600, 12*3600, 20000)  # +/- 12 h around the engineered crossing, ~4.3 s steps

    el_day = {name: np.empty_like(t_day) for name in SITES}
    for k, t in enumerate(t_day):
        raan = raan0 + orb['raan_dot']*t
        u = u0 + orb['u_dot']*t
        r_eci = eci_position(a, inc, raan, u)
        gmst = gmst0 + OMEGA_E*t
        r_ecef = eci_to_ecef(r_eci, gmst)
        for name, s in SITES.items():
            el, _, _ = look_angles(r_ecef, site_ecef[name], s['lat'], s['lon'])
            el_day[name][k] = el

    colors = {"Delhi (Alice)": "#1f77b4", "Mumbai (Bob)": "#d62728"}
    fig, ax = plt.subplots(figsize=(10, 4))
    for name in SITES:
        ax.plot(t_day/3600, el_day[name], label=name, color=colors[name], lw=1.0)
    ax.axhline(el_min, color='gray', ls='--', lw=1, label=f'min. elevation ({el_min:.0f} deg)')
    ax.set_xlabel('Time relative to engineered crossing (hours)')
    ax.set_ylabel('Elevation (deg)')
    ax.set_ylim(-5, 95)
    ax.set_title('24-hour elevation history, same fixed orbital plane -- most orbits miss the corridor')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)

    n_passes = {}
    for name in SITES:
        above = el_day[name] >= el_min
        edges = np.diff(above.astype(int))
        n_passes[name] = int((edges == 1).sum()) + (1 if above[0] else 0)
    print("Independent single-station passes above", el_min, "deg in +/-12 h window:", n_passes)
    print(f"Orbital period: {2*np.pi/state['orb']['n0']/60:.2f} min -> ~{24*60/(2*np.pi/state['orb']['n0']/60):.1f} orbits/day")
    return dict(t_day=t_day, el_day=el_day, n_passes=n_passes)


def run_full_part3_pipeline(show_plots=True):
    """Runs notebook cells 79-98 in order and returns the geometry_state dict
    that Part IV onward consumes."""
    print_sites_table()
    geodesy_ecef_check()
    df_orbit = orbit_design_table()
    print(df_orbit)
    validate_slant_range()
    state = construct_pass()
    print(state['df_pass'])
    plot_pass_geometry(state, show=show_plots)
    day_result = day_sweep(state, show=show_plots)
    state.update(day_result)
    return state
