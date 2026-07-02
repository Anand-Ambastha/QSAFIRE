# -*- coding: utf-8 -*-
"""
Figure A12 — Optimization Landscape R(eps, mu')
Figure A13 — Optimized epsilon vs Distance
Figure A14 — Optimized mu' vs Distance
Figure A15 — Optimized mu1 vs Distance
Figure A16 — Optimized mu2 vs Distance

A12 is verbatim logic from the original ``visualization.py``
(``plot_optimization_landscape``). A13-A16 expose
``SNSOptimizer.optimize_all()``'s existing outputs (never plotted
before) across a distance sweep — no optimizer logic changed.
"""

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from snstfqkd.security.decoy_analysis import DecoyStateAnalyzer
from snstfqkd.security.key_rate import SNSKeyRateCalculator
from snstfqkd.optimization.optimizer import SNSOptimizer


def figure_a12_optimization_landscape(
    channel_cls, channel_params_cls, distance: float = 500.0,
    mu1: float = 0.01, mu2: float = 0.15, misalignment: float = 0.0,
    f: float = 1.1, steps: int = 40, savepath: Optional[str] = None,
):
    """2-D heat-map of key rate over (epsilon, mu') parameter space."""
    analyzer = DecoyStateAnalyzer()
    calc = SNSKeyRateCalculator(f=f)
    ch = channel_cls(channel_params_cls(misalignment_error=misalignment))

    eps_vals = np.linspace(0.01, 0.5, steps)
    mu_vals = np.linspace(0.005, 0.4, steps)
    Z = np.zeros((steps, steps))

    s0 = ch.vacuum_yield()
    decoy = analyzer.three_intensity_decoy(
        mu1, mu2, s0,
        ch.intensity_yield(mu1, distance), ch.intensity_yield(mu2, distance),
        ch.qber_x_basis(mu1, distance),
    )

    for i, eps in enumerate(eps_vals):
        for j, mu_s in enumerate(mu_vals):
            if hasattr(ch, 'z_window_yield'):
                S_Z = ch.z_window_yield(mu_s, distance)
            else:
                S_Z = ch.intensity_yield(mu_s, distance)
            r = calc.key_rate(eps, mu_s, decoy.s1, decoy.e1_ph, S_Z, ch.qber_z_basis())
            Z[i, j] = r if r > 0 else np.nan

    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(
        np.log10(np.where(np.isnan(Z), 1e-20, Z)),
        origin="lower", aspect="auto",
        extent=[mu_vals[0], mu_vals[-1], eps_vals[0], eps_vals[-1]],
        cmap="viridis",
    )
    fig.colorbar(im, ax=ax, label=r"$\log_{10}(R)$")
    ax.set_xlabel(r"Signal intensity $\mu'$")
    ax.set_ylabel(r"Sending probability $\varepsilon$")
    ax.set_title(
        f"Key Rate Landscape at {distance} km"
        f"\n($e_a={misalignment}$, $\\mu_1={mu1}$, $\\mu_2={mu2}$)"
    )
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def _optimize_sweep(channel, distances, f, de_maxiter, de_popsize, seed):
    opt = SNSOptimizer(channel, f=f, de_maxiter=de_maxiter)
    results = []
    for d in distances:
        results.append(opt.optimize_all(
            d, method="differential_evolution",
            de_maxiter=de_maxiter, de_popsize=de_popsize, seed=seed,
        ))
    return results


def _optimized_param_figure(distances, values, ylabel, title, savepath):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(distances, values, "o-", ms=3, color="tab:purple")
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax


def figure_a13_optimized_epsilon_vs_distance(channel, distances: np.ndarray, f: float = 1.1,
                                              de_maxiter: int = 100, de_popsize: int = 10,
                                              seed: int = 42, savepath: Optional[str] = None):
    results = _optimize_sweep(channel, distances, f, de_maxiter, de_popsize, seed)
    return _optimized_param_figure(distances, [r.epsilon for r in results],
                                    r"Optimized $\varepsilon$", r"Optimized $\varepsilon$ vs Distance",
                                    savepath)


def figure_a14_optimized_mu_signal_vs_distance(channel, distances: np.ndarray, f: float = 1.1,
                                                de_maxiter: int = 100, de_popsize: int = 10,
                                                seed: int = 42, savepath: Optional[str] = None):
    results = _optimize_sweep(channel, distances, f, de_maxiter, de_popsize, seed)
    return _optimized_param_figure(distances, [r.mu_signal for r in results],
                                    r"Optimized $\mu'$", r"Optimized $\mu'$ vs Distance",
                                    savepath)


def figure_a15_optimized_mu1_vs_distance(channel, distances: np.ndarray, f: float = 1.1,
                                          de_maxiter: int = 100, de_popsize: int = 10,
                                          seed: int = 42, savepath: Optional[str] = None):
    results = _optimize_sweep(channel, distances, f, de_maxiter, de_popsize, seed)
    return _optimized_param_figure(distances, [r.mu1 for r in results],
                                    r"Optimized $\mu_1$", r"Optimized $\mu_1$ vs Distance",
                                    savepath)


def figure_a16_optimized_mu2_vs_distance(channel, distances: np.ndarray, f: float = 1.1,
                                          de_maxiter: int = 100, de_popsize: int = 10,
                                          seed: int = 42, savepath: Optional[str] = None):
    results = _optimize_sweep(channel, distances, f, de_maxiter, de_popsize, seed)
    return _optimized_param_figure(distances, [r.mu2 for r in results],
                                    r"Optimized $\mu_2$", r"Optimized $\mu_2$ vs Distance",
                                    savepath)
