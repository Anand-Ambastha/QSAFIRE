# -*- coding: utf-8 -*-
"""
Figure P7  — Accepted phase-slice fraction
Figure P8  — Rejected phase-slice fraction
Figure P9  — Effective-event fraction
Figure P10 — Protocol convergence: observable vs pulse count

P7/P8 reuse the geometric acceptance-fraction formula already exposed
in ``channel.detector.phase_slice_acceptance_fraction`` (same Eq. (1)
criterion the protocol simulator's ``sample_sliced_phases`` uses).
P9/P10 come directly from ``SNSProtocolSimulator.simulate_x_windows``'s
existing ``XWindowStats`` fields. No new physics.
"""

from typing import Optional, Sequence

import numpy as np
import matplotlib.pyplot as plt

from snstfqkd.channel.detector import phase_slice_acceptance_fraction
from snstfqkd.protocol.protocol_simulator import SNSProtocolSimulator


def figure_p7_accepted_phase_slice_fraction(lambda_slices: np.ndarray,
                                             savepath: Optional[str] = None):
    vals = [phase_slice_acceptance_fraction(l) for l in lambda_slices]
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(lambda_slices, vals, "o-", ms=3, color="tab:blue")
    ax.set_xlabel(r"Phase-slice half-width $\lambda$")
    ax.set_ylabel("Accepted fraction")
    ax.set_title("Accepted Phase-Slice Fraction")
    ax.grid(True, alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def figure_p8_rejected_phase_slice_fraction(lambda_slices: np.ndarray,
                                             savepath: Optional[str] = None):
    vals = [1.0 - phase_slice_acceptance_fraction(l) for l in lambda_slices]
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(lambda_slices, vals, "o-", ms=3, color="tab:red")
    ax.set_xlabel(r"Phase-slice half-width $\lambda$")
    ax.set_ylabel("Rejected fraction")
    ax.set_title("Rejected Phase-Slice Fraction")
    ax.grid(True, alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def figure_p9_effective_event_fraction_vs_distance(channel, distances: np.ndarray,
                                                     epsilon: float = 0.05, mu_signal: float = 0.05,
                                                     mu1: float = 0.01, mu2: float = 0.15,
                                                     n_pulses: int = 2_000_000, seed: int = 42,
                                                     savepath: Optional[str] = None):
    vals = []
    for d in distances:
        sim = SNSProtocolSimulator(channel, d, epsilon, mu_signal, mu1, mu2,
                                    n_pulses_x=n_pulses, n_pulses_z=n_pulses, seed=seed)
        stats = sim.simulate_x_windows(mu1)
        vals.append(stats.n_effective / stats.n_phase_slice if stats.n_phase_slice else 0.0)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(distances, vals, "o-", ms=3, color="tab:green")
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Effective-event fraction")
    ax.set_title("Effective-Event Fraction vs Distance")
    ax.grid(True, alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def figure_p10_protocol_convergence(channel, distance: float = 300.0,
                                     pulse_counts: Sequence[int] = (
                                         10_000, 50_000, 200_000, 1_000_000, 5_000_000, 20_000_000
                                     ),
                                     mu1: float = 0.01, epsilon: float = 0.05, mu_signal: float = 0.05,
                                     mu2: float = 0.15, seed: int = 42,
                                     savepath: Optional[str] = None):
    """S_mu1 estimate vs number of simulated pulses, converging to the
    AnalyticalChannel closed-form value (dashed reference line)."""
    vals = []
    for n in pulse_counts:
        sim = SNSProtocolSimulator(channel, distance, epsilon, mu_signal, mu1, mu2,
                                    n_pulses_x=n, n_pulses_z=n, seed=seed)
        vals.append(sim.simulate_x_windows(mu1).S_mu)

    analytic_ref = channel.intensity_yield(mu1, distance)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.semilogx(pulse_counts, vals, "o-", ms=4, color="tab:blue", label="Protocol MC estimate")
    ax.axhline(analytic_ref, color="gray", linestyle="--", label="Analytical closed-form")
    ax.set_xlabel("Number of simulated pulses")
    ax.set_ylabel(r"$S_{\mu_1}$ estimate")
    ax.set_title(f"Protocol Convergence at {distance} km")
    ax.legend()
    ax.grid(True, which="both", alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax
