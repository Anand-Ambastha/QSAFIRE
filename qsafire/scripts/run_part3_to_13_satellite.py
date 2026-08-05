#!/usr/bin/env python3
"""Reproduce Parts III-XIII of the notebook: the ground-to-satellite link
extension (Phases 2-9) -- geometry, atmosphere, turbulence, seasonal
variation, real hardware losses, coupling to the SNS-TF-QKD engine,
weather-gated annual performance, joint signal/decoy optimization, final
validation consolidation, three additional city-pair links, and the fully
asymmetric per-link protocol metrics.

Usage:
    python scripts/run_part3_to_13_satellite.py [--no-plots]
"""
import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-plots", action="store_true", help="Skip matplotlib figure rendering.")
    args = parser.parse_args()
    show_plots = not args.no_plots

    print("=" * 70 + "\nPART III -- Ground-to-Satellite Geometry Engine (Phase 2)\n" + "=" * 70)
    state = run_full_part3_pipeline(show_plots=show_plots)

    print("\n" + "=" * 70 + "\nPART IV -- Static Atmospheric Transmittance (Phase 3)\n" + "=" * 70)
    atm = run_full_part4_pipeline(state, show_plots=show_plots)

    print("\n" + "=" * 70 + "\nPART V -- Slant-Path Turbulence & Scintillation (Phase 4)\n" + "=" * 70)
    turb = run_full_part5_pipeline(state, show_plots=show_plots)

    print("\n" + "=" * 70 + "\nPART VI -- Seasonal Extension\n" + "=" * 70)
    seasonal = run_full_part6_pipeline(state, atm["tau_R_810"], show_plots=show_plots)

    print("\n" + "=" * 70 + "\nPART VII -- PAT, Optics and Detector Real-Loss Module (Phase 5)\n" + "=" * 70)
    phase5 = run_full_part7_pipeline(state, atm["atm_results"], show_plots=show_plots)

    print("\n" + "=" * 70 + "\nPART VIII -- Coupling to the SNS-TF-QKD Engine (Phase 6)\n" + "=" * 70)
    phase6 = run_full_part8_pipeline(state, atm["atm_results"], turb["turb_results_site"],
                                      phase5["phase5_results"], show_plots=show_plots)

    print("\n" + "=" * 70 + "\nPART IX -- Weather-Gated Availability & Annual Performance (Phase 7)\n" + "=" * 70)
    weather = run_full_part9_pipeline(state, atm, turb, phase5, phase6, show_plots=show_plots)

    print("\n" + "=" * 70 + "\nPART X -- Joint Signal/Decoy Optimization (Phase 8)\n" + "=" * 70)
    phase10 = run_full_part10_pipeline(state, phase5, turb, phase6, weather, show_plots=show_plots)

    print("\n" + "=" * 70 + "\nPART XI -- Final Validation Consolidation (Phase 8)\n" + "=" * 70)
    phase11 = run_full_part11_pipeline(state, atm, turb, phase5, phase6, weather, phase10, show_plots=show_plots)

    print("\n" + "=" * 70 + "\nPART XII -- Three Additional City-Pair Links (Phase 9)\n" + "=" * 70)
    links = run_full_part12_pipeline(state, show_plots=show_plots)

    print("\n" + "=" * 70 + "\nPART XIII -- Per-Link SNS-TF-QKD Metrics (Phase 9 final)\n" + "=" * 70)
    metrics = run_full_part13_pipeline(links, show_plots=show_plots)

    print("\nDone. Figures written to the working directory; see docs/workflow.md for the full "
          "cell-by-cell -> module map.")


if __name__ == "__main__":
    main()
