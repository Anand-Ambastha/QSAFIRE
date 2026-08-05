"""Part XII orchestrator: Three Additional City-Pair Links (Phase 9).

Reproduces notebook cells 225-237 in order. Consumes ``state`` (Part III, for
``orb``/``a``/``inc``/``gmst0``) and ``phase6_state`` (Part VIII, for
``eta_point_val`` indirectly via Part VII's pointing case). Returns a
``links_state`` dict.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.satellite.geometry import slant_from_elev
from src.atmosphere.transmittance import LAMBDA_UM, TAU_GAS_ZENITH, kasten_young_airmass, tau_aerosol_zenith, tau_rayleigh_zenith
from src.atmosphere.turbulence import sigma_R2_slant, cn2_hv, LAMBDA_PHASE4_M
from src.atmosphere.seasonal import A0_BASELINE
from src.pat.losses import (
    ETA_RX_INTERNAL, ETA_DET_REAL, E_A_MISALIGNMENT, R_TX_HZ, R_NOISE_BASELINE_HZ, P_DARK_REAL,
    POINTING_CASES, eta_diffraction, eta_pointing,
)
from src.optimization.joint import optimize_joint_ensemble
from src.links.city_pairs import CITY_DB, LINK_PAIRS, LINK_COLORS
from src.links.pass_builder import build_pass, evaluate_link


def print_city_db():
    """Cell 227 (script part)."""
    for name, p in CITY_DB.items():
        print(f"{name:12s}: lat={p['lat']:7.4f} lon={p['lon']:8.4f} alt={p['alt_km']*1000:6.1f} m   "
              f"AOD500={p['aod500']:.2f}  alpha={p['alpha']:.2f}  turb_mult={p['turb_mult']:.2f}   [{p['conf']}]")


def construct_all_passes(state):
    """Cell 229 (script part)."""
    orb, a, inc, gmst0 = state['orb'], state['a'], state['inc'], state['gmst0']
    pass_results = {}
    for city_a, city_b in LINK_PAIRS:
        pr = build_pass(city_a, city_b, orb, a, inc, gmst0)
        pass_results[(city_a, city_b)] = pr
        n_common = pr['common'].sum()
        print(f"{city_a:10s}-{city_b:10s}: great-circle {pr['dist_km']:7.1f} km   "
              f"common-visibility duration = {n_common*pr['dt']:6.1f} s   "
              f"max elevation (A/B) = {pr['el_a'].max():.1f}/{pr['el_b'].max():.1f} deg")
    return pass_results


def evaluate_all_links(pass_results):
    """Cell 231 (script part)."""
    eta_point_val = eta_pointing(POINTING_CASES["Conservative default ([S2025])"]['sigma_j'])
    link_summary_rows = []
    link_eval_cache = {}
    for city_a, city_b in LINK_PAIRS:
        pr = pass_results[(city_a, city_b)]
        ev = evaluate_link(city_a, city_b, pr, eta_point_val)
        link_eval_cache[(city_a, city_b)] = ev
        if ev is None:
            link_summary_rows.append(dict(link=f"{city_a}-{city_b}", distance_km=round(pr['dist_km'], 1),
                                           status="NO common-visibility window"))
            continue
        Ra, Rb = ev["A"]["R"], ev["B"]["R"]
        link_summary_rows.append(dict(
            link=f"{city_a}-{city_b}", distance_km=round(pr['dist_km'], 1),
            common_vis_s=round(pr['common'].sum()*pr['dt'], 1),
            elev_A=round(ev["A"]["elevation"], 1), elev_B=round(ev["B"]["elevation"], 1),
            eta_total_A=f"{ev['A']['eta_total']:.2e}", eta_total_B=f"{ev['B']['eta_total']:.2e}",
            R_A_kbps=round(max(Ra, 0)*R_TX_HZ/1000, 4), R_B_kbps=round(max(Rb, 0)*R_TX_HZ/1000, 4),
            status="OK" if (Ra > 0 and Rb > 0) else ("A only" if Ra > 0 else ("B only" if Rb > 0 else "NO KEY (either station)")))
        )
    df_links = pd.DataFrame(link_summary_rows)
    print(df_links)
    return df_links, link_eval_cache, eta_point_val


def plot_multilink_comparison(pass_results, link_eval_cache, save_path="phase12_multilink_comparison.png", show=True):
    """Cell 233."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for city_a, city_b in LINK_PAIRS:
        pr = pass_results[(city_a, city_b)]
        label = f"{city_a}-{city_b}"
        axes[0].plot(pr['t'], pr['el_a'], color=LINK_COLORS[label], lw=1.6, ls='-', label=f"{label} ({city_a})")
        axes[0].plot(pr['t'], pr['el_b'], color=LINK_COLORS[label], lw=1.6, ls='--', label=f"{label} ({city_b})")
    axes[0].axhline(20, color='gray', ls=':', lw=1)
    axes[0].set_xlabel('Time relative to engineered crossing (s)')
    axes[0].set_ylabel('Elevation (deg)')
    axes[0].set_title('Elevation profiles, all four links')
    axes[0].legend(fontsize=7, ncol=1)
    axes[0].grid(alpha=0.3)

    rates = []
    for city_a, city_b in LINK_PAIRS:
        ev = link_eval_cache[(city_a, city_b)]
        rates.append(max(min(ev["A"]["R"], ev["B"]["R"]), 0)*R_TX_HZ/1000 if ev else 0.0)
    labels = [f"{a}-{b}" for a, b in LINK_PAIRS]
    axes[1].bar(labels, rates, color=[LINK_COLORS[l] for l in labels])
    axes[1].set_ylabel('Bottleneck key rate (kbit/s)\n(min of the two stations)')
    axes[1].set_title('Key rate vs. link, at best-elevation point')
    axes[1].tick_params(axis='x', rotation=20)
    axes[1].grid(alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


def elevation_sweep_all_links(eta_point_val, save_path="phase12_keyrate_vs_elevation_alllinks.png", show=True):
    """Cell 235."""
    el_sweep2 = np.linspace(20, 60, 17)
    R_vs_elev_links = {}
    for city_a, city_b in LINK_PAIRS:
        label = f"{city_a}-{city_b}"
        cp_a, cp_b = CITY_DB[city_a], CITY_DB[city_b]
        curve_a, curve_b = [], []
        for el_test in el_sweep2:
            zen = 90 - el_test
            air = kasten_young_airmass(zen)
            L_m_test = slant_from_elev(el_test, 500.0)*1000.0
            for cp, curve in [(cp_a, curve_a), (cp_b, curve_b)]:
                tau_a_ = tau_aerosol_zenith(cp['aod500'], cp['alpha'], 500.0)
                tau_ray_ = tau_rayleigh_zenith(LAMBDA_UM, P_site_hpa=1013.25*np.exp(-cp['alt_km']/8.5))
                eta_atm_ = np.exp(-(tau_ray_+tau_a_+TAU_GAS_ZENITH)*air)
                eta_tot_ = eta_diffraction(L_m_test)*eta_atm_*eta_point_val*ETA_RX_INTERNAL*ETA_DET_REAL
                A0_c = A0_BASELINE*cp['turb_mult']
                sig_ = sigma_R2_slant(zen, LAMBDA_PHASE4_M, lambda h, A0=A0_c: cn2_hv(h, A0=A0))
                p_dark_ = P_DARK_REAL + (2.0*R_NOISE_BASELINE_HZ/R_TX_HZ)/1000.0
                R_, _, _, _ = optimize_joint_ensemble(eta_tot_, E_A_MISALIGNMENT, sig_, p_dark_, n_mc=6000, n_grid=7)
                curve.append(max(R_, 0)*R_TX_HZ/1000)
        R_vs_elev_links[label] = (curve_a, curve_b)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for label, (ca, cb) in R_vs_elev_links.items():
        ax.plot(el_sweep2, ca, color=LINK_COLORS[label], lw=2.0, label=f"{label.split('-')[0]} ({label})")
        ax.plot(el_sweep2, cb, color=LINK_COLORS[label], lw=1.2, ls='--')
    ax.set_xlabel('Elevation angle (deg)')
    ax.set_ylabel('Key rate (kbit/s)')
    ax.set_title('Key rate vs. elevation, all four links (solid=first city, dashed=second)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return R_vs_elev_links


def run_full_part12_pipeline(state, show_plots=True):
    """Runs notebook cells 225-237 in order. Returns links_state dict for
    Part XIII."""
    print_city_db()
    pass_results = construct_all_passes(state)
    df_links, link_eval_cache, eta_point_val = evaluate_all_links(pass_results)
    plot_multilink_comparison(pass_results, link_eval_cache, show=show_plots)
    R_vs_elev_links = elevation_sweep_all_links(eta_point_val, show=show_plots)
    return dict(pass_results=pass_results, df_links=df_links, link_eval_cache=link_eval_cache,
                eta_point_val=eta_point_val, R_vs_elev_links=R_vs_elev_links)
