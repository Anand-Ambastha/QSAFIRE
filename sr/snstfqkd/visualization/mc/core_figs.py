# -*- coding: utf-8 -*-
"""
Figure MC1 — MC Key Rate vs Distance
Figure MC2 — MC Gain vs Distance
Figure MC3 — MC S_mu vs Distance
Figure MC4 — MC QBER vs Distance
Figure MC5 — MC S_Z vs Distance

Same plotting pattern as the analytical figures, driven by a
``MonteCarloChannel`` instance instead of ``AnalyticalChannel``. No new
physics — these simply call the existing MonteCarloChannel methods.
"""

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from snstfqkd.security.decoy_analysis import DecoyStateAnalyzer
from snstfqkd.security.key_rate import SNSKeyRateCalculator
from snstfqkd.security.bounds import get_S_Z


def _curve(distances, values, ylabel, title, logy=True, color="tab:orange", savepath=None):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    values = np.asarray(values, dtype=float)
    if logy:
        ax.semilogy(distances, np.maximum(values, 1e-300), color=color, marker="o", ms=3)
    else:
        ax.plot(distances, values, color=color, marker="o", ms=3)
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, which="both" if logy else "major", alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def figure_mc1_key_rate_vs_distance(mc_channel, distances: np.ndarray,
                                     epsilon: float = 0.05, mu_signal: float = 0.05,
                                     mu1: float = 0.01, mu2: float = 0.15, f: float = 1.1,
                                     savepath: Optional[str] = None):
    calc = SNSKeyRateCalculator(f=f)
    results = calc.key_rate_vs_distance(distances, mc_channel, epsilon, mu_signal, mu1, mu2)
    rates = [r.key_rate for r in results]
    return _curve(distances, rates, "Key rate per pulse", "MC Key Rate vs Distance",
                  savepath=savepath)


def figure_mc2_gain_vs_distance(mc_channel, distances: np.ndarray, mu_signal: float = 0.05,
                                 savepath: Optional[str] = None):
    vals = [mc_channel.gain(mu_signal, d) for d in distances]
    return _curve(distances, vals, "Gain", "MC Gain vs Distance", savepath=savepath)


def figure_mc3_s_mu_vs_distance(mc_channel, distances: np.ndarray, mu1: float = 0.01,
                                 mu2: float = 0.15, savepath: Optional[str] = None):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    s_mu1 = [mc_channel.intensity_yield(mu1, d) for d in distances]
    s_mu2 = [mc_channel.intensity_yield(mu2, d) for d in distances]
    ax.semilogy(distances, np.maximum(s_mu1, 1e-300), "o-", ms=3, label=r"$S_{\mu_1}$")
    ax.semilogy(distances, np.maximum(s_mu2, 1e-300), "s-", ms=3, label=r"$S_{\mu_2}$")
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel(r"$S_\mu$")
    ax.set_title("MC $S_\\mu$ vs Distance")
    ax.legend()
    ax.grid(True, which="both", alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def figure_mc4_qber_vs_distance(mc_channel, distances: np.ndarray, mu1: float = 0.01,
                                 savepath: Optional[str] = None):
    vals = [mc_channel.qber_x_basis(mu1, d) for d in distances]
    return _curve(distances, vals, r"$E_{\mu_1}^X$", "MC QBER vs Distance", logy=False,
                  color="tab:red", savepath=savepath)


def figure_mc5_s_z_vs_distance(mc_channel, distances: np.ndarray, mu_signal: float = 0.05,
                                savepath: Optional[str] = None):
    vals = [get_S_Z(mc_channel, mu_signal, d) for d in distances]
    return _curve(distances, vals, "$S_Z$", "MC $S_Z$ vs Distance", savepath=savepath)
