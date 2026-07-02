# -*- coding: utf-8 -*-
"""
snstfqkd
==========

Research-grade SNS-TF-QKD simulation package.

Layers (per architecture audit):
    channel/       — BaseChannel, AnalyticalChannel, MonteCarloChannel
    protocol/      — SNSProtocolSimulator (independent Layer 3)
    security/      — decoy-state analysis (Eq. 44/45), key rate (Eq. 4)
    optimization/  — SNSOptimizer
    visualization/ — all Figure A/MC/P/C generators
    workflows/     — per-mode orchestration
    reports/       — summary_report.txt / distance_breakdown.csv
    config/        — default parameters and named presets
    cli.py, main.py — unified command-line interface
"""

from snstfqkd.security.decoy_analysis import DecoyStateAnalyzer, DecoyResult
from snstfqkd.security.key_rate import SNSKeyRateCalculator, KeyRateResult, binary_entropy
from snstfqkd.optimization.optimizer import SNSOptimizer, OptimizationResult

__all__ = [
    "DecoyStateAnalyzer", "DecoyResult",
    "SNSKeyRateCalculator", "KeyRateResult", "binary_entropy",
    "SNSOptimizer", "OptimizationResult",
]

__version__ = "2.0.0"
