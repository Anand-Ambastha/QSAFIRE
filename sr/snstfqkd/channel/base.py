# -*- coding: utf-8 -*-
"""
snstfqkd.channel.base
======================

Shared data classes and the abstract channel interface for the
SNS-TF-QKD channel layer.

This module contains everything that was IDENTICAL between the old
``sns_tfqkd.py`` (analytical) and ``sns_tfqkd_mc.py`` (Monte-Carlo)
files: ``ChannelParameters``, ``SimulationParameters``,
``DetectorModel``, and the (now-formalized) ``BaseChannel`` interface
that ``AnalyticalChannel`` and ``MonteCarloChannel`` both implement.

No equations were changed during this extraction — this is a pure
move of duplicated code into one place.

Reference: Wang, X.-B., Yu, Z.-W., & Hu, X.-L. (2018). Sending or not
sending: Twin-field quantum key distribution with large misalignment
error. arXiv:1805.09222v9.
"""

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Tuple

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# ============================================================================
# DATA CLASSES  (unchanged from sns_tfqkd.py / sns_tfqkd_mc.py — identical
# in both original files)
# ============================================================================

@dataclass
class ChannelParameters:
    """
    Physical parameters for the SNS-TF-QKD channel.

    TF-QKD geometry::

        Alice ──── Charlie ──── Bob
                  ←  L/2  →←  L/2  →

    Parameters
    ----------
    fiber_loss : float
        Fiber attenuation alpha [dB/km].  Default 0.2 dB/km (SMF-28 at 1550 nm).
    detector_efficiency : float
        Single-photon detector efficiency eta_d in (0, 1].  Default 0.80.
    dark_count_rate : float
        Dark count probability per pulse d.  Default 10^-11 (paper value).
    error_correction_efficiency : float
        Shannon correction factor f >= 1.  Default 1.1 (paper value).
    misalignment_error : float
        Single-photon interference misalignment e_a in [0, 0.5].
        Fraction of single-photon events that go to the wrong detector
        in the X-basis.  Default 0.0.
    visibility : float
        Optical interference visibility V in [0, 1].  Retained for API
        compatibility but NOT used in intensity_yield() or qber_x_basis().
        The Wang 2018 paper sets V = 1 and encodes imperfection entirely
        through e_a.  Default 0.99.
    """

    # Fiber
    fiber_loss: float = 0.2              # dB/km

    # Detector  (Wang et al. 2018 values)
    detector_efficiency: float = 0.8    # eta_d
    dark_count_rate: float = 1e-11      # d per pulse

    # Protocol
    error_correction_efficiency: float = 1.1   # f

    # Misalignment (X-basis only)
    misalignment_error: float = 0.0    # e_a

    # Interference visibility
    visibility: float = 0.99           # V

    def __post_init__(self) -> None:
        if not (0 < self.detector_efficiency <= 1):
            raise ValueError("detector_efficiency must be in (0, 1]")
        if self.dark_count_rate < 0:
            raise ValueError("dark_count_rate must be non-negative")
        if self.fiber_loss <= 0:
            raise ValueError("fiber_loss must be positive")
        if self.error_correction_efficiency < 1:
            raise ValueError("error_correction_efficiency must be >= 1")
        if not (0 <= self.misalignment_error <= 0.5):
            raise ValueError("misalignment_error must be in [0, 0.5]")
        if not (0 <= self.visibility <= 1):
            raise ValueError("visibility must be in [0, 1]")


@dataclass
class SimulationParameters:
    """Numerical sweep parameters (unchanged from original)."""

    min_distance: float = 0.0
    max_distance: float = 800.0
    distance_points: int = 100

    epsilon_range: Tuple[float, float] = (0.05, 0.95)
    mu_signal_range: Tuple[float, float] = (0.001, 0.5)
    mu_decoy_range: Tuple[float, float] = (0.0001, 0.4)

    epsilon_steps: int = 20
    mu_signal_steps: int = 20
    mu_decoy_steps: int = 15


# ============================================================================
# DETECTOR MODEL  (unchanged API, identical in both original files)
# ============================================================================

class DetectorModel:
    """Single-photon detector model (unchanged public API)."""

    def __init__(self, efficiency: float = 0.8, dark_count_rate: float = 1e-11):
        self.efficiency = efficiency
        self.dark_count_rate = dark_count_rate

    def click_probability(self, mean_photons: float) -> float:
        """P(click) = 1 - exp(-eta * n - d)."""
        return 1.0 - np.exp(-self.efficiency * mean_photons - self.dark_count_rate)

    def dark_count_probability(self) -> float:
        return self.dark_count_rate


# ============================================================================
# BASE CHANNEL INTERFACE
# ============================================================================

class BaseChannel(ABC):
    """
    Abstract interface implemented by both ``AnalyticalChannel`` and
    ``MonteCarloChannel``.

    This formalizes the duck-typed interface that ``decoy_analysis.py``,
    ``key_rate.py``, ``optimizer.py``, and ``visualization.py`` already
    relied on in the original codebase (confirmed by the ``FakeChannel``
    test double in the original ``test_phase23.py``). No behavior is
    changed by introducing this ABC — it only documents/enforces the
    contract that already existed implicitly.
    """

    params: ChannelParameters

    @abstractmethod
    def channel_transmittance(self, distance: float) -> float:
        """Fiber-only transmittance for one arm of length L/2."""

    @abstractmethod
    def total_transmittance(self, distance: float) -> float:
        """Total single-arm transmittance including detector efficiency."""

    @abstractmethod
    def vacuum_yield(self) -> float:
        """Vacuum counting rate at Charlie (mu = 0)."""

    @abstractmethod
    def intensity_yield(self, intensity: float, distance: float) -> float:
        """X-window single-click counting rate S_mu at Charlie's BS."""

    @abstractmethod
    def z_window_yield(self, intensity: float, distance: float) -> float:
        """Z-window single-click counting rate S_Z at Charlie."""

    @abstractmethod
    def gain(self, intensity: float, distance: float) -> float:
        """Counting rate at Charlie (equals intensity_yield for X-windows)."""

    @abstractmethod
    def qber_x_basis(self, intensity: float, distance: float) -> float:
        """X-basis quantum bit error rate E_mu^X."""

    @abstractmethod
    def qber_z_basis(self) -> float:
        """Z-basis quantum bit error rate E_Z."""
