# -*- coding: utf-8 -*-
"""
Figure P1 — Z-window (effective) fraction vs Distance
Figure P2 — X-window (effective) fraction vs Distance
Figure P3 — SS frequency (both Alice & Bob send in a Z-window)
Figure P4 — SN frequency (Alice sends, Bob doesn't)
Figure P5 — NS frequency (Bob sends, Alice doesn't)
Figure P6 — NN frequency (neither sends)

P1/P2 use the already-computed S_Z / S_mu observables from
``SNSProtocolSimulator``. P3-P6 expose the Alice/Bob send decision
frequencies, which the original simulator drew (in
``simulate_z_windows``) but never counted/reported — done here by
directly re-using ``protocol.pulse.sample_z_window_intensities``
(unchanged function, same RNG draw) and classifying its output. No new
physics: this only counts outcomes of a draw the simulator already made.
"""

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from snstfqkd.protocol.protocol_simulator import SNSProtocolSimulator
from snstfqkd.protocol.pulse import sample_z_window_intensities


def _curve(distances, values, ylabel, title, savepath=None, color="tab:blue"):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(distances, values, "o-", ms=3, color=color)
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(0, max(1.0, max(values) * 1.1) if len(values) else 1.0)
    ax.grid(True, alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def figure_p1_z_window_fraction_vs_distance(channel, distances: np.ndarray,
                                             epsilon: float = 0.05, mu_signal: float = 0.05,
                                             mu1: float = 0.01, mu2: float = 0.15,
                                             n_pulses: int = 2_000_000, seed: int = 42,
                                             savepath: Optional[str] = None):
    vals = []
    for d in distances:
        sim = SNSProtocolSimulator(channel, d, epsilon, mu_signal, mu1, mu2,
                                    n_pulses_x=n_pulses, n_pulses_z=n_pulses, seed=seed)
        vals.append(sim.simulate_z_windows().S_Z)
    return _curve(distances, vals, "Z-window effective fraction ($S_Z$)",
                  "Z-window Fraction vs Distance", savepath)


def figure_p2_x_window_fraction_vs_distance(channel, distances: np.ndarray,
                                             epsilon: float = 0.05, mu_signal: float = 0.05,
                                             mu1: float = 0.01, mu2: float = 0.15,
                                             n_pulses: int = 2_000_000, seed: int = 42,
                                             savepath: Optional[str] = None):
    vals = []
    for d in distances:
        sim = SNSProtocolSimulator(channel, d, epsilon, mu_signal, mu1, mu2,
                                    n_pulses_x=n_pulses, n_pulses_z=n_pulses, seed=seed)
        vals.append(sim.simulate_x_windows(mu1).S_mu)
    return _curve(distances, vals, r"X-window effective fraction ($S_{\mu_1}$)",
                  "X-window Fraction vs Distance", savepath)


def _send_frequencies(epsilon: float, mu_signal: float, n_pulses: int = 5_000_000, seed: int = 42):
    rng = np.random.default_rng(seed)
    mu_A, mu_B = sample_z_window_intensities(rng, n_pulses, epsilon, mu_signal)
    a_sends = mu_A > 0.0
    b_sends = mu_B > 0.0
    ss = float(np.mean(a_sends & b_sends))
    sn = float(np.mean(a_sends & ~b_sends))
    ns = float(np.mean(~a_sends & b_sends))
    nn = float(np.mean(~a_sends & ~b_sends))
    return ss, sn, ns, nn


def figure_p3_ss_frequency_vs_epsilon(epsilons: np.ndarray, mu_signal: float = 0.05,
                                       n_pulses: int = 2_000_000, seed: int = 42,
                                       savepath: Optional[str] = None):
    vals = [_send_frequencies(e, mu_signal, n_pulses, seed)[0] for e in epsilons]
    return _curve(epsilons, vals, "SS frequency", "SS Frequency vs $\\varepsilon$",
                  savepath=savepath, color="tab:green")


def figure_p4_sn_frequency_vs_epsilon(epsilons: np.ndarray, mu_signal: float = 0.05,
                                       n_pulses: int = 2_000_000, seed: int = 42,
                                       savepath: Optional[str] = None):
    vals = [_send_frequencies(e, mu_signal, n_pulses, seed)[1] for e in epsilons]
    return _curve(epsilons, vals, "SN frequency", "SN Frequency vs $\\varepsilon$",
                  savepath=savepath, color="tab:orange")


def figure_p5_ns_frequency_vs_epsilon(epsilons: np.ndarray, mu_signal: float = 0.05,
                                       n_pulses: int = 2_000_000, seed: int = 42,
                                       savepath: Optional[str] = None):
    vals = [_send_frequencies(e, mu_signal, n_pulses, seed)[2] for e in epsilons]
    return _curve(epsilons, vals, "NS frequency", "NS Frequency vs $\\varepsilon$",
                  savepath=savepath, color="tab:purple")


def figure_p6_nn_frequency_vs_epsilon(epsilons: np.ndarray, mu_signal: float = 0.05,
                                       n_pulses: int = 2_000_000, seed: int = 42,
                                       savepath: Optional[str] = None):
    vals = [_send_frequencies(e, mu_signal, n_pulses, seed)[3] for e in epsilons]
    return _curve(epsilons, vals, "NN frequency", "NN Frequency vs $\\varepsilon$",
                  savepath=savepath, color="tab:red")
