# -*- coding: utf-8 -*-
"""
Figure MC6  — Phase Slice Acceptance Fraction vs Distance
Figure MC7  — Effective Event Fraction vs Distance
Figure MC8  — Left Detector Click Rate vs Distance
Figure MC9  — Right Detector Click Rate vs Distance
Figure MC10 — Double Click Probability vs Distance
Figure MC11 — No Click Probability vs Distance

These quantities were not previously exposed by ``MonteCarloChannel``'s
public API (only S_mu, E_mu were returned). Per the requirement, the
underlying quantities are exposed (via
``channel.detector.mc_beamsplitter_diagnostics`` — a NEW function
reproducing the SAME physical model verbatim, just returning more of
what it already computes) without changing any physics formula.
"""

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from snstfqkd.channel.detector import mc_beamsplitter_diagnostics, _derive_mc_seed


def _diagnostics_sweep(mc_channel, distances, intensity):
    out = []
    for d in distances:
        seed = _derive_mc_seed(mc_channel._mc_seed, intensity, d, mc_channel.params.misalignment_error)
        diag = mc_beamsplitter_diagnostics(
            mc_channel.params.fiber_loss,
            mc_channel.params.detector_efficiency,
            mc_channel.params.dark_count_rate,
            mc_channel.params.misalignment_error,
            float(intensity),
            float(d),
            mc_channel._mc_n_samples,
            mc_channel._mc_lambda_slice,
            seed,
        )
        out.append(diag)
    return out


def _diag_curve(distances, values, ylabel, title, savepath=None):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(distances, values, "o-", ms=3, color="tab:green")
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def figure_mc6_phase_slice_acceptance_vs_distance(mc_channel, distances: np.ndarray,
                                                    intensity: float = 0.05,
                                                    savepath: Optional[str] = None):
    diags = _diagnostics_sweep(mc_channel, distances, intensity)
    vals = [d["phase_slice_acceptance"] for d in diags]
    return _diag_curve(distances, vals, "Phase-slice acceptance fraction",
                        "Phase Slice Acceptance Fraction vs Distance", savepath)


def figure_mc7_effective_event_fraction_vs_distance(mc_channel, distances: np.ndarray,
                                                      intensity: float = 0.05,
                                                      savepath: Optional[str] = None):
    diags = _diagnostics_sweep(mc_channel, distances, intensity)
    vals = [d["effective_fraction"] for d in diags]
    return _diag_curve(distances, vals, "Effective event fraction",
                        "Effective Event Fraction vs Distance", savepath)


def figure_mc8_left_click_rate_vs_distance(mc_channel, distances: np.ndarray,
                                            intensity: float = 0.05,
                                            savepath: Optional[str] = None):
    diags = _diagnostics_sweep(mc_channel, distances, intensity)
    vals = [d["P_L_mean"] for d in diags]
    return _diag_curve(distances, vals, "Left detector click rate",
                        "Left Detector Click Rate vs Distance", savepath)


def figure_mc9_right_click_rate_vs_distance(mc_channel, distances: np.ndarray,
                                             intensity: float = 0.05,
                                             savepath: Optional[str] = None):
    diags = _diagnostics_sweep(mc_channel, distances, intensity)
    vals = [d["P_R_mean"] for d in diags]
    return _diag_curve(distances, vals, "Right detector click rate",
                        "Right Detector Click Rate vs Distance", savepath)


def figure_mc10_double_click_probability_vs_distance(mc_channel, distances: np.ndarray,
                                                       intensity: float = 0.05,
                                                       savepath: Optional[str] = None):
    diags = _diagnostics_sweep(mc_channel, distances, intensity)
    vals = [d["double_click_fraction"] for d in diags]
    return _diag_curve(distances, vals, "Double-click probability",
                        "Double Click Probability vs Distance", savepath)


def figure_mc11_no_click_probability_vs_distance(mc_channel, distances: np.ndarray,
                                                   intensity: float = 0.05,
                                                   savepath: Optional[str] = None):
    diags = _diagnostics_sweep(mc_channel, distances, intensity)
    vals = [d["no_click_fraction"] for d in diags]
    return _diag_curve(distances, vals, "No-click probability",
                        "No Click Probability vs Distance", savepath)
