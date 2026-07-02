# -*- coding: utf-8 -*-
"""
snstfqkd.workflows.comparison_workflow
=========================================

Drives ``AnalyticalChannel`` vs ``MonteCarloChannel`` through Figures
C1-C8 for the ``--mode compare`` CLI target, and writes a comparison
``summary_report.txt`` / ``distance_breakdown.csv`` with both channels'
values side by side.
"""

import os
import time

import numpy as np
import pandas as pd

from snstfqkd.channel.analytical_channel import AnalyticalChannel
from snstfqkd.channel.mc_channel import MonteCarloChannel
from snstfqkd.config.presets import get_preset
from snstfqkd.visualization.comparison import (
    figure_c1_s_mu_comparison, figure_c2_qber_comparison, figure_c3_s_z_comparison,
    figure_c4_gain_comparison, figure_c5_key_rate_comparison,
    figure_c6_relative_error_s_mu, figure_c7_relative_error_qber,
    figure_c8_relative_error_key_rate,
)


def run_comparison_workflow(
    output_root: str = "results/comparison",
    preset: str = "paper",
    distances: np.ndarray = None,
    fast: bool = True,
):
    t0 = time.time()
    fig_dir = os.path.join(output_root, "figures")
    csv_dir = os.path.join(output_root, "csv")
    rep_dir = os.path.join(output_root, "reports")
    for d in (fig_dir, csv_dir, rep_dir):
        os.makedirs(d, exist_ok=True)

    params, cfg = get_preset(preset)
    analytical = AnalyticalChannel(params)
    mc = MonteCarloChannel(params, mc_n_samples=200_000 if fast else 2_000_000)

    if distances is None:
        distances = np.array(sorted(set(cfg["distances"])), dtype=float)
    sweep_distances = np.linspace(distances.min(), distances.max(), 12 if fast else 30)

    _, _, s_mu_a, s_mu_m = figure_c1_s_mu_comparison(
        analytical, mc, sweep_distances, mu1=cfg["mu1"],
        savepath=os.path.join(fig_dir, "fig_c1_s_mu_comparison.png"))
    _, _, qber_a, qber_m = figure_c2_qber_comparison(
        analytical, mc, sweep_distances, mu1=cfg["mu1"],
        savepath=os.path.join(fig_dir, "fig_c2_qber_comparison.png"))
    _, _, sz_a, sz_m = figure_c3_s_z_comparison(
        analytical, mc, sweep_distances, mu_signal=cfg["mu_signal"],
        savepath=os.path.join(fig_dir, "fig_c3_s_z_comparison.png"))
    _, _, gain_a, gain_m = figure_c4_gain_comparison(
        analytical, mc, sweep_distances, mu_signal=cfg["mu_signal"],
        savepath=os.path.join(fig_dir, "fig_c4_gain_comparison.png"))
    _, _, rate_a, rate_m = figure_c5_key_rate_comparison(
        analytical, mc, sweep_distances, epsilon=cfg["epsilon"], mu_signal=cfg["mu_signal"],
        mu1=cfg["mu1"], mu2=cfg["mu2"],
        savepath=os.path.join(fig_dir, "fig_c5_key_rate_comparison.png"))
    figure_c6_relative_error_s_mu(
        analytical, mc, sweep_distances, mu1=cfg["mu1"],
        savepath=os.path.join(fig_dir, "fig_c6_relative_error_s_mu.png"))
    figure_c7_relative_error_qber(
        analytical, mc, sweep_distances, mu1=cfg["mu1"],
        savepath=os.path.join(fig_dir, "fig_c7_relative_error_qber.png"))
    figure_c8_relative_error_key_rate(
        analytical, mc, sweep_distances, epsilon=cfg["epsilon"], mu_signal=cfg["mu_signal"],
        mu1=cfg["mu1"], mu2=cfg["mu2"],
        savepath=os.path.join(fig_dir, "fig_c8_relative_error_key_rate.png"))

    df = pd.DataFrame({
        "distance_km": sweep_distances,
        "S_mu1_analytical": s_mu_a, "S_mu1_mc": s_mu_m,
        "E_mu1_analytical": qber_a, "E_mu1_mc": qber_m,
        "S_Z_analytical": sz_a, "S_Z_mc": sz_m,
        "gain_analytical": gain_a, "gain_mc": gain_m,
        "R_analytical": rate_a, "R_mc": rate_m,
    })
    df.to_csv(os.path.join(csv_dir, "distance_breakdown.csv"), index=False)

    peak_a = float(np.max(rate_a)) if len(rate_a) else 0.0
    peak_m = float(np.max(rate_m)) if len(rate_m) else 0.0
    lines = [
        "=" * 72,
        "SNS-TF-QKD Summary Report — mode: compare",
        "=" * 72,
        "",
        f"Peak key rate (Analytical): {peak_a:.6e}",
        f"Peak key rate (Monte Carlo): {peak_m:.6e}",
        f"Max relative error in key rate: "
        f"{float(np.max(np.abs(np.array(rate_m) - np.array(rate_a)) / np.maximum(np.abs(rate_a), 1e-300))):.4%}",
        f"Runtime: {time.time() - t0:.2f} s",
        "Simulation mode: compare",
        "",
    ]
    with open(os.path.join(rep_dir, "summary_report.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))

    return df
