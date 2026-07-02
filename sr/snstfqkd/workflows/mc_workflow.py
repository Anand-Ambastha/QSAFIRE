# -*- coding: utf-8 -*-
"""
snstfqkd.workflows.mc_workflow
=================================

Drives ``MonteCarloChannel`` through Figures MC1-MC12, the
``distance_breakdown.csv`` and ``summary_report.txt`` for the
``--mode mc`` CLI target.
"""

import os
import time

import numpy as np

from snstfqkd.channel.mc_channel import MonteCarloChannel
from snstfqkd.config.presets import get_preset
from snstfqkd.reports.report_generator import (
    generate_distance_breakdown, write_breakdown_csv, write_summary_report,
)
from snstfqkd.visualization.mc import (
    figure_mc1_key_rate_vs_distance, figure_mc2_gain_vs_distance,
    figure_mc3_s_mu_vs_distance, figure_mc4_qber_vs_distance,
    figure_mc5_s_z_vs_distance, figure_mc6_phase_slice_acceptance_vs_distance,
    figure_mc7_effective_event_fraction_vs_distance, figure_mc8_left_click_rate_vs_distance,
    figure_mc9_right_click_rate_vs_distance, figure_mc10_double_click_probability_vs_distance,
    figure_mc11_no_click_probability_vs_distance, figure_mc12_runtime_vs_samples,
)


def run_mc_workflow(
    output_root: str = "results/mc",
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
    n_samples = 200_000 if fast else 2_000_000
    channel = MonteCarloChannel(params, mc_n_samples=n_samples)

    if distances is None:
        distances = np.array(sorted(set(cfg["distances"])), dtype=float)
    sweep_distances = np.linspace(distances.min(), distances.max(), 15 if fast else 40)

    figure_mc1_key_rate_vs_distance(
        channel, sweep_distances, epsilon=cfg["epsilon"], mu_signal=cfg["mu_signal"],
        mu1=cfg["mu1"], mu2=cfg["mu2"],
        savepath=os.path.join(fig_dir, "fig_mc1_key_rate_vs_distance.png"))
    figure_mc2_gain_vs_distance(
        channel, sweep_distances, mu_signal=cfg["mu_signal"],
        savepath=os.path.join(fig_dir, "fig_mc2_gain_vs_distance.png"))
    figure_mc3_s_mu_vs_distance(
        channel, sweep_distances, mu1=cfg["mu1"], mu2=cfg["mu2"],
        savepath=os.path.join(fig_dir, "fig_mc3_s_mu_vs_distance.png"))
    figure_mc4_qber_vs_distance(
        channel, sweep_distances, mu1=cfg["mu1"],
        savepath=os.path.join(fig_dir, "fig_mc4_qber_vs_distance.png"))
    figure_mc5_s_z_vs_distance(
        channel, sweep_distances, mu_signal=cfg["mu_signal"],
        savepath=os.path.join(fig_dir, "fig_mc5_s_z_vs_distance.png"))
    figure_mc6_phase_slice_acceptance_vs_distance(
        channel, sweep_distances, intensity=cfg["mu_signal"],
        savepath=os.path.join(fig_dir, "fig_mc6_phase_slice_acceptance_vs_distance.png"))
    figure_mc7_effective_event_fraction_vs_distance(
        channel, sweep_distances, intensity=cfg["mu_signal"],
        savepath=os.path.join(fig_dir, "fig_mc7_effective_event_fraction_vs_distance.png"))
    figure_mc8_left_click_rate_vs_distance(
        channel, sweep_distances, intensity=cfg["mu_signal"],
        savepath=os.path.join(fig_dir, "fig_mc8_left_click_rate_vs_distance.png"))
    figure_mc9_right_click_rate_vs_distance(
        channel, sweep_distances, intensity=cfg["mu_signal"],
        savepath=os.path.join(fig_dir, "fig_mc9_right_click_rate_vs_distance.png"))
    figure_mc10_double_click_probability_vs_distance(
        channel, sweep_distances, intensity=cfg["mu_signal"],
        savepath=os.path.join(fig_dir, "fig_mc10_double_click_probability_vs_distance.png"))
    figure_mc11_no_click_probability_vs_distance(
        channel, sweep_distances, intensity=cfg["mu_signal"],
        savepath=os.path.join(fig_dir, "fig_mc11_no_click_probability_vs_distance.png"))
    figure_mc12_runtime_vs_samples(
        channel, savepath=os.path.join(fig_dir, "fig_mc12_runtime_vs_samples.png"))

    df = generate_distance_breakdown(
        channel, distances,
        de_maxiter=40 if fast else 150, de_popsize=8, optimize=True,
    )
    write_breakdown_csv(df, os.path.join(csv_dir, "distance_breakdown.csv"))
    write_summary_report(
        df, os.path.join(rep_dir, "summary_report.txt"),
        mode="mc", channel_params=params, runtime_seconds=time.time() - t0,
    )
    return df
