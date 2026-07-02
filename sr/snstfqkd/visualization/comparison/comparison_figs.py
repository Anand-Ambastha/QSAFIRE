# -*- coding: utf-8 -*-
"""
Figure C1 — Analytical vs MC: S_mu
Figure C2 — Analytical vs MC: QBER
Figure C3 — Analytical vs MC: S_Z
Figure C4 — Analytical vs MC: Gain
Figure C5 — Analytical vs MC: Key Rate
Figure C6 — Relative Error: S_mu
Figure C7 — Relative Error: QBER
Figure C8 — Relative Error: Key Rate

Cross-validates ``AnalyticalChannel`` against ``MonteCarloChannel`` at
matched operating points — this is the cross-check the audit flagged
as missing (Risk #1). No physics changed; this only calls existing
methods on both channel implementations and plots/differences the
results.
"""

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from snstfqkd.security.key_rate import SNSKeyRateCalculator
from snstfqkd.security.bounds import get_S_Z


def _overlay(distances, analytic_vals, mc_vals, ylabel, title, logy=True, savepath=None):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    analytic_vals = np.asarray(analytic_vals, dtype=float)
    mc_vals = np.asarray(mc_vals, dtype=float)
    plot = ax.semilogy if logy else ax.plot
    if logy:
        ax.semilogy(distances, np.maximum(analytic_vals, 1e-300), "-", color="tab:blue", label="Analytical")
        ax.semilogy(distances, np.maximum(mc_vals, 1e-300), "o", ms=4, color="tab:orange", label="Monte Carlo")
    else:
        ax.plot(distances, analytic_vals, "-", color="tab:blue", label="Analytical")
        ax.plot(distances, mc_vals, "o", ms=4, color="tab:orange", label="Monte Carlo")
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, which="both" if logy else "major", alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def _relative_error(distances, analytic_vals, mc_vals, ylabel, title, savepath=None):
    analytic_vals = np.asarray(analytic_vals, dtype=float)
    mc_vals = np.asarray(mc_vals, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_err = np.where(
            np.abs(analytic_vals) > 1e-300,
            np.abs(mc_vals - analytic_vals) / np.abs(analytic_vals),
            0.0,
        )
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.semilogy(distances, np.maximum(rel_err, 1e-12), "o-", ms=4, color="tab:red")
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def figure_c1_s_mu_comparison(analytical_channel, mc_channel, distances: np.ndarray,
                               mu1: float = 0.01, savepath: Optional[str] = None):
    a = [analytical_channel.intensity_yield(mu1, d) for d in distances]
    m = [mc_channel.intensity_yield(mu1, d) for d in distances]
    fig, ax = _overlay(distances, a, m, r"$S_{\mu_1}$", r"Analytical vs MC: $S_\mu$", savepath=savepath)
    return fig, ax, a, m


def figure_c2_qber_comparison(analytical_channel, mc_channel, distances: np.ndarray,
                               mu1: float = 0.01, savepath: Optional[str] = None):
    a = [analytical_channel.qber_x_basis(mu1, d) for d in distances]
    m = [mc_channel.qber_x_basis(mu1, d) for d in distances]
    fig, ax = _overlay(distances, a, m, r"$E_{\mu_1}^X$", "Analytical vs MC: QBER", logy=False, savepath=savepath)
    return fig, ax, a, m


def figure_c3_s_z_comparison(analytical_channel, mc_channel, distances: np.ndarray,
                              mu_signal: float = 0.05, savepath: Optional[str] = None):
    a = [get_S_Z(analytical_channel, mu_signal, d) for d in distances]
    m = [get_S_Z(mc_channel, mu_signal, d) for d in distances]
    fig, ax = _overlay(distances, a, m, "$S_Z$", "Analytical vs MC: $S_Z$", savepath=savepath)
    return fig, ax, a, m


def figure_c4_gain_comparison(analytical_channel, mc_channel, distances: np.ndarray,
                               mu_signal: float = 0.05, savepath: Optional[str] = None):
    a = [analytical_channel.gain(mu_signal, d) for d in distances]
    m = [mc_channel.gain(mu_signal, d) for d in distances]
    fig, ax = _overlay(distances, a, m, "Gain", "Analytical vs MC: Gain", savepath=savepath)
    return fig, ax, a, m


def figure_c5_key_rate_comparison(analytical_channel, mc_channel, distances: np.ndarray,
                                   epsilon: float = 0.05, mu_signal: float = 0.05,
                                   mu1: float = 0.01, mu2: float = 0.15, f: float = 1.1,
                                   savepath: Optional[str] = None):
    calc = SNSKeyRateCalculator(f=f)
    a = [r.key_rate for r in calc.key_rate_vs_distance(distances, analytical_channel, epsilon, mu_signal, mu1, mu2)]
    m = [r.key_rate for r in calc.key_rate_vs_distance(distances, mc_channel, epsilon, mu_signal, mu1, mu2)]
    fig, ax = _overlay(distances, a, m, "Key rate per pulse", "Analytical vs MC: Key Rate", savepath=savepath)
    return fig, ax, a, m


def figure_c6_relative_error_s_mu(analytical_channel, mc_channel, distances: np.ndarray,
                                   mu1: float = 0.01, savepath: Optional[str] = None):
    a = [analytical_channel.intensity_yield(mu1, d) for d in distances]
    m = [mc_channel.intensity_yield(mu1, d) for d in distances]
    return _relative_error(distances, a, m, r"Relative error in $S_{\mu_1}$",
                            "Relative Error: $S_\\mu$", savepath=savepath)


def figure_c7_relative_error_qber(analytical_channel, mc_channel, distances: np.ndarray,
                                   mu1: float = 0.01, savepath: Optional[str] = None):
    a = [analytical_channel.qber_x_basis(mu1, d) for d in distances]
    m = [mc_channel.qber_x_basis(mu1, d) for d in distances]
    return _relative_error(distances, a, m, "Relative error in QBER",
                            "Relative Error: QBER", savepath=savepath)


def figure_c8_relative_error_key_rate(analytical_channel, mc_channel, distances: np.ndarray,
                                       epsilon: float = 0.05, mu_signal: float = 0.05,
                                       mu1: float = 0.01, mu2: float = 0.15, f: float = 1.1,
                                       savepath: Optional[str] = None):
    calc = SNSKeyRateCalculator(f=f)
    a = [r.key_rate for r in calc.key_rate_vs_distance(distances, analytical_channel, epsilon, mu_signal, mu1, mu2)]
    m = [r.key_rate for r in calc.key_rate_vs_distance(distances, mc_channel, epsilon, mu_signal, mu1, mu2)]
    return _relative_error(distances, a, m, "Relative error in key rate",
                            "Relative Error: Key Rate", savepath=savepath)
