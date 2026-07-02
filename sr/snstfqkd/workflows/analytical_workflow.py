# -*- coding: utf-8 -*-
"""
snstfqkd.workflows.analytical_workflow
=========================================

Drives ``AnalyticalChannel`` through all Figure A1-A16, the
``distance_breakdown.csv``, and ``summary_report.txt`` for the
``--mode analytical`` CLI target. Pure orchestration — no physics.
"""

import os
import time

import numpy as np

from snstfqkd.channel.analytical_channel import AnalyticalChannel
from snstfqkd.config.presets import get_preset
from snstfqkd.config import defaults as D
from snstfqkd.reports.report_generator import (
    generate_distance_breakdown, write_breakdown_csv, write_summary_report,
)
from snstfqkd.visualization.analytical import (
    figure_a1_key_rate_vs_distance, figure_a2_key_rate_vs_misalignment,
    figure_a3_s1_vs_distance, figure_a4_e1ph_vs_distance,
    figure_a5_gain_vs_distance, figure_a6_transmittance_vs_distance,
    figure_a7_vacuum_yield_vs_distance, figure_a8_s_mu1_vs_distance,
    figure_a9_s_mu2_vs_distance, figure_a10_e_mu1_vs_distance,
    figure_a11_s_z_vs_distance, figure_a12_optimization_landscape,
    figure_a13_optimized_epsilon_vs_distance, figure_a14_optimized_mu_signal_vs_distance,
    figure_a15_optimized_mu1_vs_distance, figure_a16_optimized_mu2_vs_distance,
)


def run_analytical_workflow(
    output_root: str = "results/analytical",
    preset: str = "paper",
    distances: np.ndarray = None,
    fast: bool = True,
):
    """
    Generate all analytical outputs (figures, CSV, report).

    fast : bool
        If True, uses reduced optimizer iterations / sweep resolution
        so the full CLI run completes quickly. Set False for
        publication-quality settings (matches original breakdown.py
        defaults).
    """
    t0 = time.time()
    fig_dir = os.path.join(output_root, "figures")
    csv_dir = os.path.join(output_root, "csv")
    rep_dir = os.path.join(output_root, "reports")
    for d in (fig_dir, csv_dir, rep_dir):
        os.makedirs(d, exist_ok=True)

    params, cfg = get_preset(preset)
    channel = AnalyticalChannel(params)

    if distances is None:
        distances = np.array(sorted(set(cfg["distances"])), dtype=float)
    sweep_distances = np.linspace(distances.min(), distances.max(), 60 if fast else 200)

    de_maxiter = 60 if fast else D.DEFAULT_DE_MAXITER
    de_popsize = 8 if fast else D.DEFAULT_DE_POPSIZE

    figure_a1_key_rate_vs_distance(
        AnalyticalChannel, type(params), sweep_distances,
        mu1=cfg["mu1"], mu2=cfg["mu2"], epsilon=cfg["epsilon"], mu_signal=cfg["mu_signal"],
        savepath=os.path.join(fig_dir, "fig_a1_key_rate_vs_distance.png"),
        optimize=not fast, de_maxiter=de_maxiter, de_popsize=de_popsize,
    )
    figure_a2_key_rate_vs_misalignment(
        AnalyticalChannel, type(params), np.linspace(0.0, 0.45, 20),
        distance=500.0, mu1=cfg["mu1"], mu2=cfg["mu2"],
        epsilon=cfg["epsilon"], mu_signal=cfg["mu_signal"],
        savepath=os.path.join(fig_dir, "fig_a2_key_rate_vs_misalignment.png"),
        optimize=not fast, de_maxiter=de_maxiter, de_popsize=de_popsize,
    )
    figure_a3_s1_vs_distance(
        AnalyticalChannel, type(params), sweep_distances, mu1=cfg["mu1"], mu2=cfg["mu2"],
        savepath=os.path.join(fig_dir, "fig_a3_s1_vs_distance.png"),
    )
    figure_a4_e1ph_vs_distance(
        AnalyticalChannel, type(params), sweep_distances, mu1=cfg["mu1"], mu2=cfg["mu2"],
        savepath=os.path.join(fig_dir, "fig_a4_e1ph_vs_distance.png"),
    )
    figure_a5_gain_vs_distance(channel, sweep_distances, intensity=cfg["mu_signal"],
                                savepath=os.path.join(fig_dir, "fig_a5_gain_vs_distance.png"))
    figure_a6_transmittance_vs_distance(channel, sweep_distances,
                                         savepath=os.path.join(fig_dir, "fig_a6_transmittance_vs_distance.png"))
    figure_a7_vacuum_yield_vs_distance(channel, sweep_distances,
                                        savepath=os.path.join(fig_dir, "fig_a7_vacuum_yield_vs_distance.png"))
    figure_a8_s_mu1_vs_distance(channel, sweep_distances, mu1=cfg["mu1"],
                                 savepath=os.path.join(fig_dir, "fig_a8_s_mu1_vs_distance.png"))
    figure_a9_s_mu2_vs_distance(channel, sweep_distances, mu2=cfg["mu2"],
                                 savepath=os.path.join(fig_dir, "fig_a9_s_mu2_vs_distance.png"))
    figure_a10_e_mu1_vs_distance(channel, sweep_distances, mu1=cfg["mu1"],
                                  savepath=os.path.join(fig_dir, "fig_a10_e_mu1_vs_distance.png"))
    figure_a11_s_z_vs_distance(channel, sweep_distances, mu_signal=cfg["mu_signal"],
                                savepath=os.path.join(fig_dir, "fig_a11_s_z_vs_distance.png"))
    figure_a12_optimization_landscape(
        AnalyticalChannel, type(params), distance=500.0, mu1=cfg["mu1"], mu2=cfg["mu2"],
        steps=25 if fast else 40,
        savepath=os.path.join(fig_dir, "fig_a12_optimization_landscape.png"),
    )
    opt_distances = np.linspace(distances.min(), distances.max(), 12 if fast else 30)
    figure_a13_optimized_epsilon_vs_distance(
        channel, opt_distances, de_maxiter=de_maxiter, de_popsize=de_popsize,
        savepath=os.path.join(fig_dir, "fig_a13_optimized_epsilon_vs_distance.png"))
    figure_a14_optimized_mu_signal_vs_distance(
        channel, opt_distances, de_maxiter=de_maxiter, de_popsize=de_popsize,
        savepath=os.path.join(fig_dir, "fig_a14_optimized_mu_signal_vs_distance.png"))
    figure_a15_optimized_mu1_vs_distance(
        channel, opt_distances, de_maxiter=de_maxiter, de_popsize=de_popsize,
        savepath=os.path.join(fig_dir, "fig_a15_optimized_mu1_vs_distance.png"))
    figure_a16_optimized_mu2_vs_distance(
        channel, opt_distances, de_maxiter=de_maxiter, de_popsize=de_popsize,
        savepath=os.path.join(fig_dir, "fig_a16_optimized_mu2_vs_distance.png"))

    df = generate_distance_breakdown(
        channel, distances, de_maxiter=de_maxiter, de_popsize=de_popsize, optimize=True,
    )
    write_breakdown_csv(df, os.path.join(csv_dir, "distance_breakdown.csv"))
    write_summary_report(
        df, os.path.join(rep_dir, "summary_report.txt"),
        mode="analytical", channel_params=params, runtime_seconds=time.time() - t0,
    )
    return df
