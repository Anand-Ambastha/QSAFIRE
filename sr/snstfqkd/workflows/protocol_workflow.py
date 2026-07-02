# -*- coding: utf-8 -*-
"""
snstfqkd.workflows.protocol_workflow
=======================================

Drives ``SNSProtocolSimulator`` (backed by an ``AnalyticalChannel`` for
its physical parameters/transmittance, per the audit — the protocol
simulator was never coupled to ``MonteCarloChannel`` in the original
codebase) through Figures P1-P10 for the ``--mode protocol`` CLI
target.
"""

import os
import time

import numpy as np

from snstfqkd.channel.analytical_channel import AnalyticalChannel
from snstfqkd.config.presets import get_preset
from snstfqkd.protocol.protocol_simulator import SNSProtocolSimulator, evaluate_key_rate_from_simulation
from snstfqkd.security.decoy_analysis import DecoyStateAnalyzer
from snstfqkd.security.key_rate import SNSKeyRateCalculator
from snstfqkd.reports.report_generator import write_summary_report
import pandas as pd

from snstfqkd.visualization.protocol import (
    figure_p1_z_window_fraction_vs_distance, figure_p2_x_window_fraction_vs_distance,
    figure_p3_ss_frequency_vs_epsilon, figure_p4_sn_frequency_vs_epsilon,
    figure_p5_ns_frequency_vs_epsilon, figure_p6_nn_frequency_vs_epsilon,
    figure_p7_accepted_phase_slice_fraction, figure_p8_rejected_phase_slice_fraction,
    figure_p9_effective_event_fraction_vs_distance, figure_p10_protocol_convergence,
)


def run_protocol_workflow(
    output_root: str = "results/protocol",
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
    channel = AnalyticalChannel(params)

    if distances is None:
        distances = np.array(sorted(set(cfg["distances"])), dtype=float)
    sweep_distances = np.linspace(distances.min(), min(distances.max(), 600.0), 8 if fast else 15)
    n_pulses = 500_000 if fast else 5_000_000

    figure_p1_z_window_fraction_vs_distance(
        channel, sweep_distances, epsilon=cfg["epsilon"], mu_signal=cfg["mu_signal"],
        mu1=cfg["mu1"], mu2=cfg["mu2"], n_pulses=n_pulses,
        savepath=os.path.join(fig_dir, "fig_p1_z_window_fraction.png"))
    figure_p2_x_window_fraction_vs_distance(
        channel, sweep_distances, epsilon=cfg["epsilon"], mu_signal=cfg["mu_signal"],
        mu1=cfg["mu1"], mu2=cfg["mu2"], n_pulses=n_pulses,
        savepath=os.path.join(fig_dir, "fig_p2_x_window_fraction.png"))

    epsilons = np.linspace(0.01, 0.5, 15)
    figure_p3_ss_frequency_vs_epsilon(epsilons, mu_signal=cfg["mu_signal"], n_pulses=n_pulses,
                                       savepath=os.path.join(fig_dir, "fig_p3_ss_frequency.png"))
    figure_p4_sn_frequency_vs_epsilon(epsilons, mu_signal=cfg["mu_signal"], n_pulses=n_pulses,
                                       savepath=os.path.join(fig_dir, "fig_p4_sn_frequency.png"))
    figure_p5_ns_frequency_vs_epsilon(epsilons, mu_signal=cfg["mu_signal"], n_pulses=n_pulses,
                                       savepath=os.path.join(fig_dir, "fig_p5_ns_frequency.png"))
    figure_p6_nn_frequency_vs_epsilon(epsilons, mu_signal=cfg["mu_signal"], n_pulses=n_pulses,
                                       savepath=os.path.join(fig_dir, "fig_p6_nn_frequency.png"))

    lambda_slices = np.logspace(-5, -1, 15)
    figure_p7_accepted_phase_slice_fraction(
        lambda_slices, savepath=os.path.join(fig_dir, "fig_p7_accepted_phase_slice_fraction.png"))
    figure_p8_rejected_phase_slice_fraction(
        lambda_slices, savepath=os.path.join(fig_dir, "fig_p8_rejected_phase_slice_fraction.png"))
    figure_p9_effective_event_fraction_vs_distance(
        channel, sweep_distances, epsilon=cfg["epsilon"], mu_signal=cfg["mu_signal"],
        mu1=cfg["mu1"], mu2=cfg["mu2"], n_pulses=n_pulses,
        savepath=os.path.join(fig_dir, "fig_p9_effective_event_fraction.png"))
    figure_p10_protocol_convergence(
        channel, distance=300.0, mu1=cfg["mu1"], epsilon=cfg["epsilon"],
        mu_signal=cfg["mu_signal"], mu2=cfg["mu2"],
        savepath=os.path.join(fig_dir, "fig_p10_protocol_convergence.png"))

    # distance_breakdown.csv from protocol-simulated observables
    analyzer = DecoyStateAnalyzer()
    calc = SNSKeyRateCalculator()
    rows = []
    for d in distances:
        sim = SNSProtocolSimulator(
            channel, d, epsilon=cfg["epsilon"], mu_signal=cfg["mu_signal"],
            mu1=cfg["mu1"], mu2=cfg["mu2"], n_pulses_x=n_pulses, n_pulses_z=n_pulses,
        )
        stats = sim.run()
        decoy, rate = evaluate_key_rate_from_simulation(stats, analyzer, calc)
        rows.append(dict(
            distance_km=d, epsilon=cfg["epsilon"], mu_signal=cfg["mu_signal"],
            mu1=cfg["mu1"], mu2=cfg["mu2"], S_mu0=stats.S_mu0, S_mu1=stats.S_mu1,
            S_mu2=stats.S_mu2, E_mu1=stats.E_X_mu1, s1=decoy.s1, e1ph=decoy.e1_ph,
            S_Z=stats.S_Z, R=rate,
        ))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(csv_dir, "distance_breakdown.csv"), index=False)

    write_summary_report(
        df.rename(columns={"R": "R"}).assign(
            gain=df["S_mu1"], positive_term=df["R"], correction_term=0.0,
        ) if len(df) else df,
        os.path.join(rep_dir, "summary_report.txt"),
        mode="protocol", channel_params=params, runtime_seconds=time.time() - t0,
    )
    return df
