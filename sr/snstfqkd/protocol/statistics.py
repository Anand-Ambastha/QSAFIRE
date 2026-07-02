# -*- coding: utf-8 -*-
"""
snstfqkd.protocol.statistics
==============================

Result containers for the protocol-level Monte Carlo simulator
(``SNSProtocolSimulator``). Verbatim move from the original
``protocol_simulator.py`` — no fields, defaults, or semantics changed.
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class XWindowStats:
    """Observed X-window statistics for one intensity setting."""
    intensity: float
    n_trials: int
    n_phase_slice: int          # number of pulse pairs passing Eq. (1).
                                 # Equal to n_trials by construction, since
                                 # phases are now drawn DIRECTLY from the
                                 # Eq. (1)-accepted region (see
                                 # sample_sliced_phases) rather than via
                                 # rejection sampling. Retained as a
                                 # separate field for API stability and to
                                 # keep S_mu's definition (n_effective /
                                 # n_phase_slice) explicit and self-
                                 # documenting, matching the paper's
                                 # "effective-event rate of X_mu-windows".
    n_effective: int            # single-click effective events
    n_errors: int                # wrong-detector clicks among effective events
    n_double_click: int
    n_no_click: int
    S_mu: float                  # observed yield = n_effective / n_phase_slice
    E_mu: float                  # observed QBER = n_errors / n_effective


@dataclass
class ZWindowStats:
    """Observed Z-window statistics."""
    mu_signal: float
    n_trials: int
    n_z_windows: int             # both committed signal window
    n_effective: int             # single-click effective events
    S_Z: float                   # observed yield = n_effective / n_z_windows


@dataclass
class ProtocolStatistics:
    """
    Full set of Monte-Carlo-observed protocol statistics, in the exact
    format consumed (unmodified) by DecoyStateAnalyzer and
    SNSKeyRateCalculator.
    """
    distance: float
    mu1: float
    mu2: float
    mu_signal: float
    epsilon: float
    lambda_slice: float
    seed: int

    S_mu0: float      # vacuum yield (X-window, intensity 0)
    S_mu1: float      # X-window yield at mu1            -> Eq. (44) input
    S_mu2: float      # X-window yield at mu2            -> Eq. (44) input
    E_X_mu1: float    # X-window QBER at mu1             -> Eq. (45) input
    S_Z: float        # Z-window yield at mu_signal      -> Eq. (4)  input
    E_Z: float = 0.0  # Z-basis QBER (= 0 in SNS, Wang Sec. III)

    x_window_details: dict = field(default_factory=dict)   # mu -> XWindowStats
    z_window_details: Optional[ZWindowStats] = None
