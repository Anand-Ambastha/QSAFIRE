"""Part VI orchestrator: Explicit Winter/Summer/Monsoon Seasonal Extension.

Reproduces notebook cells 131-142 in order. Consumes ``state`` (Part III),
``tau_R_810``/``SITE_TAU_ZENITH`` inputs (Part IV), and
``sigma_R2_slant``/``cn2_hv`` (Part V). Returns a ``seasonal_state`` dict
carrying ``seasonal_results`` forward.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.satellite.geometry import SITES
from src.atmosphere.transmittance import tau_aerosol_zenith, kasten_young_airmass, TAU_GAS_ZENITH
from src.atmosphere.turbulence import cn2_hv, sigma_R2_slant, LAMBDA_PHASE4_M
from src.atmosphere.seasonal import SEASONAL_AEROSOL, SEASONAL_TURBULENCE, A0_BASELINE, SEASONS
from src.fso.channel_models import gamma_gamma_params, classify_turbulence_regime

COLORS = {"Delhi (Alice)": "#1f77b4", "Mumbai (Bob)": "#d62728"}


def seasonal_aerosol_table():
    """Cell 135 (script part)."""
    rows = []
    for station, seasons in SEASONAL_AEROSOL.items():
        for season, p in seasons.items():
            tau_a = tau_aerosol_zenith(p['aod0'], p['alpha'], p['lam0_nm'])
            rows.append(dict(station=station, season=season, aod=p['aod0'], alpha=p['alpha'],
                              tau_aerosol_810nm=round(tau_a, 3), confidence=p['conf']))
    return pd.DataFrame(rows)


def seasonal_turbulence_table():
    """Cell 137 (script part)."""
    from src.atmosphere.turbulence import fried_r0
    rows = []
    for station, seasons in SEASONAL_TURBULENCE.items():
        for season, p in seasons.items():
            r0z, _ = fried_r0(0.0, LAMBDA_PHASE4_M, A0=p['A0'])
            rows.append(dict(station=station, season=season, A0=f"{p['A0']:.2e}",
                              zenith_r0_cm=round(r0z*100, 2), confidence=p['conf']))
    return pd.DataFrame(rows)


def recompute_pass_seasonal(state, tau_R_810):
    """Cell 139: recompute the Phase 2 pass geometry season by season,
    holding the orbital geometry fixed."""
    results, idx = state['results'], state['idx']

    seasonal_results = {}
    for station in SITES:
        zenith_deg_full = 90.0 - results[station]['el']
        airmass_full = kasten_young_airmass(zenith_deg_full)
        for season in SEASONS:
            aero = SEASONAL_AEROSOL[station][season]
            turb = SEASONAL_TURBULENCE[station][season]

            tau_a = tau_aerosol_zenith(aero['aod0'], aero['alpha'], aero['lam0_nm'])
            tau_zenith_total = tau_R_810 + tau_a + TAU_GAS_ZENITH
            tau_slant = tau_zenith_total * airmass_full
            eta_atm = np.exp(-tau_slant)
            loss_dB = -10*np.log10(eta_atm)

            sigma_R2 = np.array([sigma_R2_slant(z, LAMBDA_PHASE4_M, lambda h, A0=turb['A0']: cn2_hv(h, A0=A0))
                                  for z in zenith_deg_full])
            ag, bg = gamma_gamma_params(sigma_R2)
            regime = np.asarray(classify_turbulence_regime(sigma_R2))

            seasonal_results[(station, season)] = dict(
                eta_atm=eta_atm, loss_dB=loss_dB, sigma_R2=sigma_R2, alpha_g=ag, beta_g=bg, regime=regime)

    rows = []
    for (station, season), d in seasonal_results.items():
        loss_win = d['loss_dB'][idx]
        s2_win = d['sigma_R2'][idx]
        regimes, counts = np.unique(d['regime'][idx], return_counts=True)
        regime_str = ", ".join(f"{r}:{100*c/len(idx):.0f}%" for r, c in zip(regimes, counts))
        rows.append(dict(station=station, season=season,
                          atm_loss_dB_range=f"[{loss_win.min():.2f}, {loss_win.max():.2f}]",
                          sigma_R2_range=f"[{s2_win.min():.3f}, {s2_win.max():.3f}]",
                          turbulence_regime=regime_str))
    df_seasonal_summary = pd.DataFrame(rows)
    return seasonal_results, df_seasonal_summary


def plot_seasonal_comparison(state, seasonal_results, save_path="phase5_seasonal_comparison.png", show=True):
    """Cell 140."""
    t_grid, t_start, t_end = state['t_grid'], state['t_start'], state['t_end']
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    season_ls = {"Winter (DJF)": "-", "Summer/pre-monsoon (MAMJ)": "--", "Monsoon (JAS)": ":"}

    for station in SITES:
        ax = axes[0, 0] if station == "Delhi (Alice)" else axes[0, 1]
        for season, ls in season_ls.items():
            ax.plot(t_grid, seasonal_results[(station, season)]['loss_dB'], color=COLORS[station],
                    ls=ls, lw=1.8, label=season)
        ax.axvspan(t_start, t_end, color='green', alpha=0.10)
        ax.set_title(f'{station}: atmospheric loss by season')
        ax.set_ylabel('Loss (dB)')
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)

    for station in SITES:
        ax = axes[1, 0] if station == "Delhi (Alice)" else axes[1, 1]
        for season, ls in season_ls.items():
            ax.plot(t_grid, seasonal_results[(station, season)]['sigma_R2'], color=COLORS[station],
                    ls=ls, lw=1.8, label=season)
        ax.axhline(1.0, color='gray', ls=':', lw=1)
        ax.axvspan(t_start, t_end, color='green', alpha=0.10)
        ax.set_title(f'{station}: Rytov variance by season')
        ax.set_ylabel(r'$\sigma_R^2$')
        ax.set_xlabel('Time relative to engineered crossing (s)')
        ax.set_ylim(0, 3)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


def run_full_part6_pipeline(state, tau_R_810, show_plots=True):
    """Runs notebook cells 131-142 in order. Returns seasonal_state dict."""
    df_aero = seasonal_aerosol_table()
    print(df_aero)
    df_turb = seasonal_turbulence_table()
    print(df_turb)
    seasonal_results, df_seasonal_summary = recompute_pass_seasonal(state, tau_R_810)
    print(df_seasonal_summary)
    plot_seasonal_comparison(state, seasonal_results, show=show_plots)
    return dict(seasonal_results=seasonal_results, df_seasonal_summary=df_seasonal_summary)
