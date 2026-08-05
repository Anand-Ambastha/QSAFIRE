"""Interactive Streamlit dashboard for the SNS-TF-QKD satellite-link pipeline
(Parts III-XIII: geometry -> atmosphere -> turbulence -> real hardware losses
-> key rate -> weather-gated availability -> joint optimization -> multi-city
links -> final per-link metrics).

Run with:
    streamlit run dashboard/app.py

The full pipeline involves Monte Carlo ensemble optimization at many points
and takes on the order of a few minutes on first load; the result is cached
for the rest of the session (see ``load_pipeline()`` below).
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.satellite.geometry import SITES
from src.links.city_pairs import CITY_DB, LINK_PAIRS
from src.experiments.satellite_experiments import run_full_part3_pipeline
from src.experiments.atmosphere_experiments import run_full_part4_pipeline
from src.experiments.turbulence_experiments import run_full_part5_pipeline
from src.experiments.seasonal_experiments import run_full_part6_pipeline
from src.experiments.pat_experiments import run_full_part7_pipeline
from src.experiments.phase6_experiments import run_full_part8_pipeline
from src.experiments.weather_experiments import run_full_part9_pipeline
from src.experiments.phase10_experiments import run_full_part10_pipeline
from src.experiments.consolidation_experiments import run_full_part11_pipeline
from src.experiments.city_links_experiments import run_full_part12_pipeline
from src.experiments.link_metrics_experiments import run_full_part13_pipeline

STATION_COLORS = {"Delhi (Alice)": "#1f77b4", "Mumbai (Bob)": "#d62728"}

st.set_page_config(page_title="QSAFire — SNS-TF-QKD Satellite Link", layout="wide")


# --------------------------------------------------------------------------
# Pipeline execution (cached — computed once per session)
# --------------------------------------------------------------------------

@st.cache_data(show_spinner=False)
def load_pipeline():
    """Runs the full Parts III-XIII pipeline (no plotting) and returns every
    stage's result dict, keyed by stage name."""
    state = run_full_part3_pipeline(show_plots=False)
    atm = run_full_part4_pipeline(state, show_plots=False)
    turb = run_full_part5_pipeline(state, show_plots=False)
    seasonal = run_full_part6_pipeline(state, atm["tau_R_810"], show_plots=False)
    phase5 = run_full_part7_pipeline(state, atm["atm_results"], show_plots=False)
    phase6 = run_full_part8_pipeline(state, atm["atm_results"], turb["turb_results_site"],
                                      phase5["phase5_results"], show_plots=False)
    weather = run_full_part9_pipeline(state, atm, turb, phase5, phase6, show_plots=False)
    phase10 = run_full_part10_pipeline(state, phase5, turb, phase6, weather, show_plots=False)
    phase11 = run_full_part11_pipeline(state, atm, turb, phase5, phase6, weather, phase10, show_plots=False)
    links = run_full_part12_pipeline(state, show_plots=False)
    metrics = run_full_part13_pipeline(links, show_plots=False)
    return dict(state=state, atm=atm, turb=turb, seasonal=seasonal, phase5=phase5,
                phase6=phase6, weather=weather, phase10=phase10, phase11=phase11,
                links=links, metrics=metrics)


def new_fig(figsize=(9, 4.5)):
    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax


# --------------------------------------------------------------------------
# Sidebar
# --------------------------------------------------------------------------

st.sidebar.title("QSAFire")
st.sidebar.caption("Satellite-Assisted SNS-TF-QKD — Delhi-Mumbai and three "
                    "additional Indian city-pair links.")
st.sidebar.markdown(
    "Ground-to-satellite orbital geometry, atmospheric transmittance, "
    "slant-path turbulence, real hardware/background-light losses, "
    "weather-gated annual availability, and joint signal/decoy "
    "optimization, feeding a fully asymmetric SNS-TF-QKD key-rate "
    "calculation."
)
run_clicked = st.sidebar.button("▶ Run / refresh full simulation", type="primary")
if run_clicked:
    load_pipeline.clear()
st.sidebar.caption("First run takes a few minutes (Monte Carlo ensemble "
                    "optimization); cached for the rest of the session.")

with st.spinner("Running the full satellite-link pipeline (geometry → atmosphere → "
                 "turbulence → real losses → key rate → weather → optimization → "
                 "links → final metrics)... this can take a few minutes on first load."):
    data = load_pipeline()

state = data["state"]
atm = data["atm"]
turb = data["turb"]
phase5 = data["phase5"]
phase6 = data["phase6"]
weather = data["weather"]
phase10 = data["phase10"]
phase11 = data["phase11"]
links = data["links"]
metrics = data["metrics"]

t_grid, results, idx = state["t_grid"], state["results"], state["idx"]
t_start, t_end = state["t_start"], state["t_end"]

st.title("QSAFire — Satellite-Assisted SNS-TF-QKD Dashboard")

# --------------------------------------------------------------------------
# Headline metrics
# --------------------------------------------------------------------------

dm_metrics = metrics["link_metrics_asym"][("Delhi", "Mumbai")]
dm_skr_kbps = max(dm_metrics["R"], 0) * 100e6 / 1000  # R_TX_HZ = 100e6, see src/pat/losses.py
best_link = max(metrics["link_metrics_asym"].items(), key=lambda kv: max(kv[1]["R"], 0))
best_link_name = f"{best_link[0][0]}-{best_link[0][1]}"
best_link_skr = max(best_link[1]["R"], 0) * 100e6 / 1000

c1, c2, c3, c4 = st.columns(4)
c1.metric("Delhi-Mumbai distance", f"{links['pass_results'][('Delhi', 'Mumbai')]['dist_km']:.1f} km")
c2.metric("Delhi-Mumbai SKR", f"{dm_skr_kbps:.4f} kbit/s")
c3.metric("Best link", best_link_name, f"{best_link_skr:.4f} kbit/s")
c4.metric("Common-visibility window", f"{state['common'].sum()*(t_grid[1]-t_grid[0]):.1f} s")

st.divider()

tabs = st.tabs([
    "Pass Geometry", "Atmosphere", "Turbulence", "Real-Loss Key Rate",
    "Weather & Annual", "Multi-City Links", "Final Metrics",
])

# --------------------------------------------------------------------------
# Tab 1: Pass geometry
# --------------------------------------------------------------------------

with tabs[0]:
    st.subheader("Delhi-Mumbai dual-downlink pass")
    st.caption("A 500 km sun-synchronous satellite pass engineered to cross the "
               "Delhi-Mumbai great-circle midpoint, with both stations simultaneously "
               "above a 20° minimum elevation.")

    col1, col2 = st.columns(2)
    with col1:
        fig, ax = new_fig()
        for name in SITES:
            ax.plot(t_grid, results[name]["el"], label=name, color=STATION_COLORS[name], lw=1.8)
        ax.axhline(state["el_min"], color="gray", ls="--", lw=1, label="min. elevation")
        ax.axvspan(t_start, t_end, color="green", alpha=0.12, label="common-visibility window")
        ax.set_xlabel("Time relative to crossing (s)"); ax.set_ylabel("Elevation (deg)")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        st.pyplot(fig); plt.close(fig)
    with col2:
        fig, ax = new_fig()
        for name in SITES:
            ax.plot(t_grid, results[name]["rng"], label=name, color=STATION_COLORS[name], lw=1.8)
        ax.axvspan(t_start, t_end, color="green", alpha=0.12)
        ax.set_xlabel("Time relative to crossing (s)"); ax.set_ylabel("Slant range (km)")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        st.pyplot(fig); plt.close(fig)

    st.dataframe(state["df_pass"], use_container_width=True)

    st.subheader("Ground sites")
    site_rows = [dict(city=k.split(" (")[0], role=k, lat=v["lat"], lon=v["lon"], alt_km=v["alt_km"])
                 for k, v in SITES.items()]
    st.map(pd.DataFrame(site_rows).rename(columns={"lat": "latitude", "lon": "longitude"}),
           latitude="latitude", longitude="longitude", size=20000)

# --------------------------------------------------------------------------
# Tab 2: Atmosphere
# --------------------------------------------------------------------------

with tabs[1]:
    st.subheader("Static atmospheric transmittance")
    st.caption("Rayleigh + aerosol + gaseous zenith optical depth, scaled to the "
               "slant path via the Kasten & Young airmass formula.")

    zenith_df = pd.DataFrame({name: {"zenith optical depth": tau,
                                      "zenith transmittance (%)": 100*np.exp(-tau)}
                               for name, tau in atm["SITE_TAU_ZENITH"].items()}).T
    st.dataframe(zenith_df, use_container_width=True)

    fig, ax = new_fig()
    for name in SITES:
        ax.plot(t_grid, atm["atm_results"][name]["eta_atm"]*100, label=name, color=STATION_COLORS[name], lw=1.8)
    ax.axvspan(t_start, t_end, color="green", alpha=0.12)
    ax.set_xlabel("Time relative to crossing (s)"); ax.set_ylabel("Atmospheric transmittance (%)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    st.pyplot(fig); plt.close(fig)

# --------------------------------------------------------------------------
# Tab 3: Turbulence
# --------------------------------------------------------------------------

with tabs[2]:
    st.subheader("Slant-path turbulence (Rytov variance, site case)")
    st.caption("Hufnagel-Valley turbulence profile, corrected spherical-wave "
               "(uplink) Rytov variance. Shared literature baseline vs. a "
               "Delhi urban-heat-island sensitivity case.")

    fig, ax = new_fig()
    for name in SITES:
        ax.plot(t_grid, turb["turb_results"][name]["sigma_R2"], color=STATION_COLORS[name],
                lw=1.2, ls="--", label=f"{name} (shared baseline)")
        ax.plot(t_grid, turb["turb_results_site"][name]["sigma_R2"], color=STATION_COLORS[name],
                lw=2.0, label=f"{name} (site case)")
    ax.axhline(1.0, color="gray", ls=":", lw=1)
    ax.axvspan(t_start, t_end, color="green", alpha=0.12)
    ax.set_xlabel("Time relative to crossing (s)"); ax.set_ylabel(r"Rytov variance $\sigma_R^2$")
    ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3)
    st.pyplot(fig); plt.close(fig)

    rows = []
    for name in SITES:
        for label, store in [("shared baseline", turb["gg_results"]), ("site case", turb["gg_results_site"])]:
            regimes, counts = np.unique(store[name]["regime"][idx], return_counts=True)
            rows.append(dict(station=name, case=label,
                              regime_mix=", ".join(f"{r}:{100*c/len(idx):.0f}%" for r, c in zip(regimes, counts)),
                              mean_alpha_g=round(store[name]["alpha_g"][idx].mean(), 2),
                              mean_beta_g=round(store[name]["beta_g"][idx].mean(), 2)))
    st.dataframe(pd.DataFrame(rows), use_container_width=True)

# --------------------------------------------------------------------------
# Tab 4: Real-loss key rate
# --------------------------------------------------------------------------

with tabs[3]:
    st.subheader("Key rate over the pass: fixed-decoy vs. joint optimization")
    st.caption("Comparison of the ensemble-optimized key rate (Part VIII, fixed "
               "20% decoy ratio) against the joint signal+decoy-intensity "
               "optimizer (Part X).")

    fig, ax = new_fig()
    for name in SITES:
        ax.plot(phase6["phase6_results"][name]["t"],
                np.clip(phase6["phase6_results"][name]["R"], 0, None)*100e6/1000,
                color=STATION_COLORS[name], lw=1.3, ls="--", label=f"{name} (fixed decoy)")
        ax.plot(phase10["phase10_results"][name]["t"],
                np.clip(phase10["phase10_results"][name]["R"], 0, None)*100e6/1000,
                color=STATION_COLORS[name], lw=2.0, label=f"{name} (joint opt.)")
    ax.set_xlabel("Time relative to crossing (s)"); ax.set_ylabel("Key rate (kbit/s)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    st.pyplot(fig); plt.close(fig)

    st.subheader("Root-cause: dark counts vs. night-sky background")
    st.caption("At each station's best-elevation point (misalignment e_a=0): "
               "dark-count-only channels stay open; adding realistic night "
               "background alone closes them.")
    st.dataframe(data["metrics"]["df_asym"][["link", "eta_A", "eta_B", "QBER_Ez_pct",
                                              "phase_err_e1ph_pct", "SKR_kbps"]],
                 use_container_width=True)

    st.subheader("Key rate vs. elevation angle (corrected)")
    el_sweep = np.linspace(20, 60, len(next(iter(phase11["R_vs_elev"].values()))))
    fig, ax = new_fig()
    for name in SITES:
        ax.plot(el_sweep, phase11["R_vs_elev"][name], color=STATION_COLORS[name], lw=2.0,
                marker="o", ms=3, label=name)
    ax.set_xlabel("Elevation angle (deg)"); ax.set_ylabel("Key rate (kbit/s)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    st.pyplot(fig); plt.close(fig)

# --------------------------------------------------------------------------
# Tab 5: Weather & annual
# --------------------------------------------------------------------------

with tabs[4]:
    st.subheader("Annual expected key volume (weather-gated, corrected)")
    st.caption("Real IMD rainy-day + winter-fog climatology, gating the orbital "
               "pass frequency; season-scaled with the joint-optimized rate.")
    st.dataframe(phase11["df_annual_v2"], use_container_width=True)

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Seasonal robustness: old vs. corrected optimizer")
        rows = []
        for station in SITES:
            for season in ["Winter", "Summer", "Monsoon"]:
                rows.append(dict(station=station, season=season,
                                  old_scale=round(weather["season_scale"].get((station, season), 0.0), 3),
                                  corrected_scale=round(phase11["season_scale_v2"].get((station, season), 0.0), 3)))
        st.dataframe(pd.DataFrame(rows), use_container_width=True)
    with col2:
        st.subheader("Sensitivity analysis (Mumbai, best-elevation point)")
        st.dataframe(weather["df_sens"], use_container_width=True)

# --------------------------------------------------------------------------
# Tab 6: Multi-city links
# --------------------------------------------------------------------------

with tabs[5]:
    st.subheader("Four candidate city-pair links")
    st.dataframe(links["df_links"], use_container_width=True)

    fig, ax = new_fig()
    labels = [f"{a}-{b}" for a, b in LINK_PAIRS]
    rates = []
    for city_a, city_b in LINK_PAIRS:
        ev = links["link_eval_cache"][(city_a, city_b)]
        rates.append(max(min(ev["A"]["R"], ev["B"]["R"]), 0)*100e6/1000 if ev else 0.0)
    bar_colors = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]
    ax.bar(labels, rates, color=bar_colors)
    ax.set_ylabel("Bottleneck key rate (kbit/s)"); ax.tick_params(axis="x", rotation=20)
    ax.grid(alpha=0.3, axis="y")
    st.pyplot(fig); plt.close(fig)

    st.subheader("Ground sites")
    city_rows = [dict(city=k, lat=v["lat"], lon=v["lon"]) for k, v in CITY_DB.items()]
    st.map(pd.DataFrame(city_rows).rename(columns={"lat": "latitude", "lon": "longitude"}),
           latitude="latitude", longitude="longitude", size=15000)

# --------------------------------------------------------------------------
# Tab 7: Final metrics
# --------------------------------------------------------------------------

with tabs[6]:
    st.subheader("Per-link SNS-TF-QKD metrics (fully asymmetric protocol)")
    st.caption("Each arm's transmittance and turbulence kept separate through "
               "gain/QBER/phase-error; validated to reduce exactly to the "
               "symmetric protocol when eta_A = eta_B.")
    st.dataframe(metrics["df_asym"], use_container_width=True)

    col1, col2, col3 = st.columns(3)
    labels = [f"{a}-{b}" for a, b in LINK_PAIRS]
    bar_colors = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]
    with col1:
        fig, ax = new_fig((5, 4))
        ax.bar(labels, [metrics["link_metrics_asym"][p]["Ez"]*100 for p in LINK_PAIRS], color=bar_colors)
        ax.set_ylabel("QBER (%)"); ax.tick_params(axis="x", rotation=30)
        st.pyplot(fig); plt.close(fig)
    with col2:
        fig, ax = new_fig((5, 4))
        ax.bar(labels, [metrics["link_metrics_asym"][p]["e1ph"]*100 for p in LINK_PAIRS], color=bar_colors)
        ax.set_ylabel("Phase-error rate (%)"); ax.tick_params(axis="x", rotation=30)
        st.pyplot(fig); plt.close(fig)
    with col3:
        fig, ax = new_fig((5, 4))
        ax.bar(labels, [max(metrics["link_metrics_asym"][p]["R"], 0)*100e6/1000 for p in LINK_PAIRS], color=bar_colors)
        ax.set_ylabel("SKR (kbit/s)"); ax.tick_params(axis="x", rotation=30)
        st.pyplot(fig); plt.close(fig)

st.divider()
st.caption("QSAFire — Quantum Satellite Analysis Framework for Integrated Research "
           "and Evaluation. See docs/theory.md and docs/api.md for the physics and "
           "function reference behind every figure on this page.")
