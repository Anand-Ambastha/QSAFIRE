# snstfqkd/security/decoy_analysis.py
# Decoy-state analysis for SNS-TF-QKD.
# Equations: Wang, Yu, Hu (2018) arXiv:1805.09222 — Eqs. 6, 8, 9, 44, 45.
# Verbatim move from the original top-level decoy_analysis.py — no
# equation changed.
import math
import numpy as np
from dataclasses import dataclass
from typing import Optional


@dataclass
class DecoyResult:
    s1: float
    e1_ph: float
    s0: float
    valid: bool
    msg: str = ""


class DecoyStateAnalyzer:
    """
    Two-mode photon-number probabilities and decoy-state bounds
    for the SNS protocol.
    """

    def p0(self, mu: float) -> float:
        """Vacuum probability. SNS Eq. (6): p0 = exp(-2μ)."""
        return np.exp(-2.0 * mu)

    def p1(self, mu: float) -> float:
        """Single-photon probability. SNS Eq. (8): p1 = 2μ exp(-2μ)."""
        return 2.0 * mu * np.exp(-2.0 * mu)

    def p2(self, mu: float) -> float:
        """Two-photon probability. SNS Eq. (9): p2 = 2μ² exp(-2μ)."""
        return 2.0 * mu ** 2 * np.exp(-2.0 * mu)

    def pk(self, mu: float, k: int) -> float:
        """General two-mode photon-number probability: e^{-2μ}(2μ)^k / k!"""
        if k < 0:
            return 0.0
        return np.exp(-2.0 * mu) * (2.0 * mu) ** k / float(math.factorial(k))

    def estimate_s1(
        self,
        mu1: float,
        mu2: float,
        S_mu1: float,
        S_mu2: float,
        s0: float,
    ) -> float:
        """Lower bound on single-photon yield. SNS Eq. (44)."""
        p2_1, p2_2 = self.p2(mu1), self.p2(mu2)
        p1_1, p1_2 = self.p1(mu1), self.p1(mu2)
        p0_1, p0_2 = self.p0(mu1), self.p0(mu2)

        num = p2_2 * (S_mu1 - p0_1 * s0) - p2_1 * (S_mu2 - p0_2 * s0)
        den = p2_2 * p1_1 - p2_1 * p1_2
        if abs(den) < 1e-20:
            return 0.0
        return num / den

    def estimate_phase_error(
        self,
        mu1: float,
        S_mu1: float,
        E_mu1: float,
        s0: float,
        s1: float,
    ) -> float:
        """Upper bound on phase-flip error rate. SNS Eq. (45)."""
        if s1 <= 0:
            return 0.5
        num = S_mu1 * E_mu1 - np.exp(-2.0 * mu1) * s0 / 2.0
        den = 2.0 * mu1 * np.exp(-2.0 * mu1) * s1
        if abs(den) < 1e-30:
            return 0.5
        return num / den

    def three_intensity_decoy(
        self,
        mu1: float,
        mu2: float,
        S_mu0: float,
        S_mu1: float,
        S_mu2: float,
        E_mu1: float,
    ) -> DecoyResult:
        """Full three-intensity decoy analysis (vacuum + mu1 + mu2)."""
        s0 = S_mu0
        s1 = float(np.clip(self.estimate_s1(mu1, mu2, S_mu1, S_mu2, s0), 0.0, 1.0))
        e1_ph = float(np.clip(
            self.estimate_phase_error(mu1, S_mu1, E_mu1, s0, s1), 0.0, 0.5
        ))
        valid = s1 > 0 and e1_ph < 0.5
        msg = "" if valid else (
            "s1 <= 0" if s1 <= 0 else "e1_ph >= 0.5"
        )
        return DecoyResult(s1=s1, e1_ph=e1_ph, s0=s0, valid=valid, msg=msg)

    def validate_bounds(self, result: DecoyResult) -> bool:
        return 0.0 <= result.s1 <= 1.0 and 0.0 <= result.e1_ph <= 0.5