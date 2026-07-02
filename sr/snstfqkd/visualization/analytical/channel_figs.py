# -*- coding: utf-8 -*-
"""
Figure A5  — Gain vs Distance
Figure A6  — Channel Transmittance vs Distance
Figure A7  — Vacuum Yield vs Distance
Figure A8  — S_mu1 vs Distance
Figure A9  — S_mu2 vs Distance
Figure A10 — E_mu1 vs Distance
Figure A11 — S_Z vs Distance

These channel observables (``gain``, ``total_transmittance``,
``vacuum_yield``, ``intensity_yield``, ``qber_x_basis``,
``z_window_yield``) already exist on every ``BaseChannel``
implementation — this module only adds plotting support for them, as
required. No underlying formula was added or changed.
"""

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt


def _single_curve_figure(distances, values, ylabel, title, logy=True, color="tab:blue",
                          savepath: Optional[str] = None):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    values = np.asarray(values, dtype=float)
    if logy:
        ax.semilogy(distances, np.maximum(values, 1e-300), color=color)
    else:
        ax.plot(distances, values, color=color)
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, which="both" if logy else "major", alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def figure_a5_gain_vs_distance(channel, distances: np.ndarray, intensity: float = 0.05,
                                savepath: Optional[str] = None):
    """Gain (channel.gain) vs distance."""
    vals = [channel.gain(intensity, d) for d in distances]
    return _single_curve_figure(distances, vals, "Gain", "Gain vs Distance", savepath=savepath)


def figure_a6_transmittance_vs_distance(channel, distances: np.ndarray,
                                         savepath: Optional[str] = None):
    """Total single-arm transmittance vs distance."""
    vals = [channel.total_transmittance(d) for d in distances]
    return _single_curve_figure(distances, vals, "Transmittance",
                                 "Channel Transmittance vs Distance", savepath=savepath)


def figure_a7_vacuum_yield_vs_distance(channel, distances: np.ndarray,
                                        savepath: Optional[str] = None):
    """Vacuum yield is distance-independent (mu=0); plotted as a flat reference line."""
    y0 = channel.vacuum_yield()
    vals = [y0 for _ in distances]
    return _single_curve_figure(distances, vals, "Vacuum yield $Y_0$",
                                 "Vacuum Yield vs Distance", savepath=savepath)


def figure_a8_s_mu1_vs_distance(channel, distances: np.ndarray, mu1: float = 0.01,
                                 savepath: Optional[str] = None):
    """X-window yield S_mu1 vs distance."""
    vals = [channel.intensity_yield(mu1, d) for d in distances]
    return _single_curve_figure(distances, vals, r"$S_{\mu_1}$", r"$S_{\mu_1}$ vs Distance",
                                 savepath=savepath)


def figure_a9_s_mu2_vs_distance(channel, distances: np.ndarray, mu2: float = 0.15,
                                 savepath: Optional[str] = None):
    """X-window yield S_mu2 vs distance."""
    vals = [channel.intensity_yield(mu2, d) for d in distances]
    return _single_curve_figure(distances, vals, r"$S_{\mu_2}$", r"$S_{\mu_2}$ vs Distance",
                                 savepath=savepath)


def figure_a10_e_mu1_vs_distance(channel, distances: np.ndarray, mu1: float = 0.01,
                                  savepath: Optional[str] = None):
    """X-basis QBER E_mu1 vs distance."""
    vals = [channel.qber_x_basis(mu1, d) for d in distances]
    return _single_curve_figure(distances, vals, r"$E_{\mu_1}^X$", r"$E_{\mu_1}^X$ vs Distance",
                                 logy=False, color="tab:red", savepath=savepath)


def figure_a11_s_z_vs_distance(channel, distances: np.ndarray, mu_signal: float = 0.05,
                                savepath: Optional[str] = None):
    """Z-window yield S_Z vs distance."""
    vals = [channel.z_window_yield(mu_signal, d) for d in distances]
    return _single_curve_figure(distances, vals, "$S_Z$", "$S_Z$ vs Distance", savepath=savepath)
