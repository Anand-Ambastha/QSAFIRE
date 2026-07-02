# -*- coding: utf-8 -*-
"""
snstfqkd.reports.report_generator
====================================

Generalizes the original top-level ``breakdown.py`` (which was
analytical-only, fixed to one preset) into a channel-agnostic CSV/report
generator usable by every workflow (analytical, mc, protocol, compare).

The per-row breakdown formula (positive_term / correction_term / R) is
copied verbatim from ``breakdown.py`` — it is algebraically identical
to ``SNSKeyRateCalculator.key_rate()`` (Eq. 4), just expanded inline so
the CSV can show the intermediate terms, exactly as the original did.
No physics changed.
"""

import time
from dataclasses import dataclass
from typing import Iterable, List, Optional

import numpy as np
import pandas as pd

from snstfqkd.security.decoy_analysis import DecoyStateAnalyzer
from snstfqkd.security.key_rate import SNSKeyRateCalculator, binary_entropy
from snstfqkd.security.bounds import get_S_Z
from snstfqkd.optimization.optimizer import SNSOptimizer


@dataclass
class BreakdownRow:
    distance_km: float
    epsilon: float
    mu_signal: float
    mu1: float
    mu2: float
    gain: float
    transmittance: float
    vacuum_yield: float
    s0: float
    S_mu1: float
    S_mu2: float
    E_mu1: float
    s1: float
    e1ph: float
    S_Z: float
    positive_term: float
    correction_term: float
    R: float


def generate_distance_breakdown(
    channel,
    distances: Iterable[float],
    f: float = 1.1,
    de_maxiter: int = 200,
    de_popsize: int = 12,
    seed: int = 42,
    optimize: bool = True,
    fixed_params: Optional[dict] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Build the distance_breakdown.csv table for one channel instance
    (AnalyticalChannel or MonteCarloChannel — anything implementing
    BaseChannel).

    Parameters
    ----------
    channel : BaseChannel
    distances : iterable of float
    optimize : bool
        If True (default), per-distance optimizes (epsilon, mu_signal,
        mu1, mu2) via SNSOptimizer.optimize_all() — same as the
        original breakdown.py. If False, uses fixed_params for all
        distances (fixed_params must supply epsilon, mu_signal, mu1, mu2).
    """
    analyzer = DecoyStateAnalyzer()
    calc = SNSKeyRateCalculator(f=f)
    optimizer = SNSOptimizer(channel, f=f, de_maxiter=de_maxiter) if optimize else None

    rows: List[BreakdownRow] = []

    for L in distances:
        if optimize:
            opt = optimizer.optimize_all(
                L,
                method="differential_evolution",
                de_maxiter=de_maxiter,
                de_popsize=de_popsize,
                seed=seed,
            )
            eps, mu_signal, mu1, mu2 = opt.epsilon, opt.mu_signal, opt.mu1, opt.mu2
        else:
            eps = fixed_params["epsilon"]
            mu_signal = fixed_params["mu_signal"]
            mu1 = fixed_params["mu1"]
            mu2 = fixed_params["mu2"]

        if verbose:
            print(f"distance = {L} km")

        s0 = channel.vacuum_yield()
        S_mu1 = channel.intensity_yield(mu1, L)
        S_mu2 = channel.intensity_yield(mu2, L)
        E_mu1 = channel.qber_x_basis(mu1, L)

        decoy = analyzer.three_intensity_decoy(mu1, mu2, s0, S_mu1, S_mu2, E_mu1)
        s1 = decoy.s1
        e1ph = decoy.e1_ph

        E_Z = channel.qber_z_basis()
        S_Z = get_S_Z(channel, mu_signal, L)

        positive = (
            2.0 * eps * (1.0 - eps)
            * mu_signal * np.exp(-mu_signal)
            * s1 * (1.0 - binary_entropy(e1ph))
        )
        correction = S_Z * calc.f * binary_entropy(E_Z)
        R = max(0.0, positive - correction)

        rows.append(BreakdownRow(
            distance_km=L,
            epsilon=eps,
            mu_signal=mu_signal,
            mu1=mu1,
            mu2=mu2,
            gain=channel.gain(mu_signal, L),
            transmittance=channel.total_transmittance(L),
            vacuum_yield=s0,
            s0=s0,
            S_mu1=S_mu1,
            S_mu2=S_mu2,
            E_mu1=E_mu1,
            s1=s1,
            e1ph=e1ph,
            S_Z=S_Z,
            positive_term=positive,
            correction_term=correction,
            R=R,
        ))

    return pd.DataFrame([r.__dict__ for r in rows])


def write_breakdown_csv(df: pd.DataFrame, path: str) -> None:
    df.to_csv(path, index=False)


def write_summary_report(
    df: pd.DataFrame,
    path: str,
    mode: str,
    channel_params,
    runtime_seconds: float,
) -> None:
    """
    Write summary_report.txt: parameters used, optimized parameters
    (at peak key rate), distance limit, runtime, peak key rate, peak
    gain, peak S_mu, peak S_Z, and simulation mode.
    """
    if len(df) == 0:
        peak_row = None
    else:
        peak_row = df.loc[df["R"].idxmax()]

    positive_rates = df[df["R"] > 0]
    max_distance = float(positive_rates["distance_km"].max()) if len(positive_rates) else float("nan")

    lines = []
    lines.append("=" * 72)
    lines.append(f"SNS-TF-QKD Summary Report — mode: {mode}")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Parameters used:")
    lines.append(f"  fiber_loss                  = {channel_params.fiber_loss} dB/km")
    lines.append(f"  detector_efficiency         = {channel_params.detector_efficiency}")
    lines.append(f"  dark_count_rate              = {channel_params.dark_count_rate}")
    lines.append(f"  error_correction_efficiency = {channel_params.error_correction_efficiency}")
    lines.append(f"  misalignment_error           = {channel_params.misalignment_error}")
    lines.append(f"  visibility                   = {channel_params.visibility}")
    lines.append("")
    if peak_row is not None:
        lines.append("Optimized parameters (at peak key rate):")
        lines.append(f"  distance   = {peak_row['distance_km']} km")
        lines.append(f"  epsilon    = {peak_row['epsilon']:.6g}")
        lines.append(f"  mu_signal  = {peak_row['mu_signal']:.6g}")
        lines.append(f"  mu1        = {peak_row['mu1']:.6g}")
        lines.append(f"  mu2        = {peak_row['mu2']:.6g}")
        lines.append("")
        lines.append(f"Distance limit (last positive key rate): {max_distance} km")
        lines.append(f"Peak key rate : {peak_row['R']:.6e}")
        lines.append(f"Peak gain     : {df['gain'].max():.6e}")
        lines.append(f"Peak S_mu1    : {df['S_mu1'].max():.6e}")
        lines.append(f"Peak S_Z      : {df['S_Z'].max():.6e}")
    else:
        lines.append("No data points computed.")
    lines.append("")
    lines.append(f"Runtime       : {runtime_seconds:.2f} s")
    lines.append(f"Simulation mode: {mode}")
    lines.append("")

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
