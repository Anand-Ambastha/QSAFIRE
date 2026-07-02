# -*- coding: utf-8 -*-
"""
Figure A3 — s1 vs Distance
Figure A4 — e1ph vs Distance

Verbatim logic from the original ``visualization.py``
(``plot_s1_vs_distance`` / ``plot_e1ph_vs_distance``).
"""

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from snstfqkd.security.decoy_analysis import DecoyStateAnalyzer


def figure_a3_s1_vs_distance(
    channel_cls, channel_params_cls, distances: np.ndarray,
    mu1: float = 0.01, mu2: float = 0.15, misalignment: float = 0.0,
    f: float = 1.1, savepath: Optional[str] = None,
):
    """Single-photon yield lower bound s1 vs distance."""
    analyzer = DecoyStateAnalyzer()
    ch = channel_cls(channel_params_cls(misalignment_error=misalignment))

    s1_vals = []
    for d in distances:
        s0 = ch.vacuum_yield()
        decoy = analyzer.three_intensity_decoy(
            mu1, mu2, s0,
            ch.intensity_yield(mu1, d), ch.intensity_yield(mu2, d),
            ch.qber_x_basis(mu1, d),
        )
        s1_vals.append(decoy.s1 if decoy.valid else 0.0)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.semilogy(distances, np.maximum(s1_vals, 1e-20), "g-",
                label=rf"$s_1$ ($e_a={misalignment}$)")
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel(r"$s_1$ (single-photon yield lower bound)")
    ax.set_title(r"$s_1$ vs Distance")
    ax.legend()
    ax.grid(True, which="both", alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def figure_a4_e1ph_vs_distance(
    channel_cls, channel_params_cls, distances: np.ndarray,
    mu1: float = 0.01, mu2: float = 0.15, misalignment: float = 0.0,
    f: float = 1.1, savepath: Optional[str] = None,
):
    """Phase-flip error upper bound e1_ph vs distance."""
    analyzer = DecoyStateAnalyzer()
    ch = channel_cls(channel_params_cls(misalignment_error=misalignment))

    e1_vals = []
    for d in distances:
        s0 = ch.vacuum_yield()
        decoy = analyzer.three_intensity_decoy(
            mu1, mu2, s0,
            ch.intensity_yield(mu1, d), ch.intensity_yield(mu2, d),
            ch.qber_x_basis(mu1, d),
        )
        e1_vals.append(decoy.e1_ph)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(distances, e1_vals, "r-", label=rf"$e_1^{{ph}}$ ($e_a={misalignment}$)")
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=1, label="Max = 0.5")
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel(r"Phase-flip error rate $e_1^{ph}$")
    ax.set_title(r"$e_1^{ph}$ vs Distance")
    ax.set_ylim(0, 0.55)
    ax.legend()
    ax.grid(True, alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax
