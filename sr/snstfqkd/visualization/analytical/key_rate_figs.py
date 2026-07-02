# -*- coding: utf-8 -*-
"""
Figure A1 — Key Rate vs Distance
Figure A2 — Key Rate vs Misalignment

Verbatim logic from the original top-level ``visualization.py``
(``plot_key_rate_vs_distance`` / ``plot_key_rate_vs_misalignment``), only
renamed and relocated. No formula changed.
"""

from typing import List, Optional

import numpy as np
import matplotlib.pyplot as plt

from snstfqkd.visualization._style import MISALIGNMENT_COLORS, LINESTYLES
from snstfqkd.security.key_rate import SNSKeyRateCalculator
from snstfqkd.optimization.optimizer import SNSOptimizer


def figure_a1_key_rate_vs_distance(
    channel_cls,
    channel_params_cls,
    distances: np.ndarray,
    misalignment_errors: List[float] = (0.00, 0.15, 0.25, 0.35, 0.45),
    mu1: float = 0.01,
    mu2: float = 0.15,
    epsilon: float = 0.05,
    mu_signal: float = 0.05,
    f: float = 1.1,
    savepath: Optional[str] = None,
    optimize: bool = True,
    de_maxiter: int = 200,
    de_popsize: int = 10,
    seed: int = 42,
):
    """Key rate vs distance for multiple misalignment values (paper Fig. 1)."""
    fig, ax = plt.subplots(figsize=(7, 5))

    for i, ea in enumerate(misalignment_errors):
        ch = channel_cls(channel_params_cls(misalignment_error=ea))

        if optimize:
            opt = SNSOptimizer(ch, f=f, de_maxiter=de_maxiter)
            rates = []
            for d in distances:
                res = opt.optimize_all(
                    d, method="differential_evolution",
                    de_maxiter=de_maxiter, de_popsize=de_popsize, seed=seed,
                )
                rates.append(res.key_rate)
            rates = np.array(rates)
        else:
            calc = SNSKeyRateCalculator(f=f)
            results = calc.key_rate_vs_distance(distances, ch, epsilon, mu_signal, mu1, mu2)
            rates = np.array([r.key_rate for r in results])

        mask = rates > 0
        if mask.any():
            ax.semilogy(
                distances[mask], rates[mask],
                linestyle=LINESTYLES[i % len(LINESTYLES)],
                color=MISALIGNMENT_COLORS[i % len(MISALIGNMENT_COLORS)],
                label=rf"$e_a={ea:.2f}$",
            )

    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Key rate per pulse")
    ax.set_title(
        "SNS-TF-QKD: Key Rate vs Distance\n"
        "(Wang et al. 2018, Fig. 1 analogue"
        + (", per-distance optimised)" if optimize else ", fixed params)")
    )
    ax.legend(loc="upper right")
    ax.set_xlim(distances[0], distances[-1])
    ax.grid(True, which="both", alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def figure_a2_key_rate_vs_misalignment(
    channel_cls,
    channel_params_cls,
    misalignments: np.ndarray,
    distance: float = 500.0,
    mu1: float = 0.01,
    mu2: float = 0.15,
    epsilon: float = 0.05,
    mu_signal: float = 0.05,
    f: float = 1.1,
    savepath: Optional[str] = None,
    optimize: bool = True,
    de_maxiter: int = 200,
    de_popsize: int = 10,
    seed: int = 42,
):
    """Key rate vs misalignment at fixed distance (paper Fig. 2)."""
    rates = []
    for ea in misalignments:
        ch = channel_cls(channel_params_cls(misalignment_error=ea))
        if optimize:
            opt = SNSOptimizer(ch, f=f, de_maxiter=de_maxiter)
            res = opt.optimize_all(
                distance, method="differential_evolution",
                de_maxiter=de_maxiter, de_popsize=de_popsize, seed=seed,
            )
            rates.append(res.key_rate)
        else:
            calc = SNSKeyRateCalculator(f=f)
            res = calc.optimized_key_rate(distance, ch, epsilon, mu_signal, mu1, mu2)
            rates.append(res.key_rate)
    rates = np.array(rates)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    mask = rates > 0
    ax.semilogy(misalignments[mask], rates[mask], "b-o", ms=4,
                label=f"Distance = {distance} km")
    ax.set_xlabel(r"Misalignment error $e_a$")
    ax.set_ylabel("Key rate per pulse")
    ax.set_title(
        "SNS-TF-QKD: Key Rate vs Misalignment\n"
        "(Wang et al. 2018, Fig. 2 analogue"
        + (", per-ea optimised)" if optimize else ", fixed params)")
    )
    ax.legend()
    ax.grid(True, which="both", alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax
