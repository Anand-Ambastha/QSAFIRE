# snstfqkd/security/key_rate.py
# Key rate calculation for SNS-TF-QKD.
# Equation reference: Wang, Yu, Hu (2018) arXiv:1805.09222 — Eq. (4).
#
# Verbatim move from the original top-level key_rate.py. The only
# structural change: the duplicated _get_S_Z() helper (identical in
# the original key_rate.py and optimizer.py) now lives once in
# security/bounds.py and is imported here — no formula changed.
#
# PATCH NOTES (carried over from original):
#   key_rate_vs_distance(): S_Z obtained from channel.z_window_yield()
#   rather than channel.intensity_yield().
#
#   Physical justification:
#     S_Z in Eq. (4) is the *Z-window* counting rate.  In a Z-window only
#     one party sends (intensity mu'), the other sends vacuum.  The BS
#     statistics are fundamentally different from X-windows (where both
#     sides send and interference occurs).  Using intensity_yield for S_Z
#     overestimates it by ~2x at short distance, inflating the
#     error-correction penalty S_Z f H(E_Z) and suppressing the key rate.
#
#   Reference: Wang 2018, Eq. (4); Step 3 definition of Z-window.
import numpy as np
from dataclasses import dataclass
from typing import List

from snstfqkd.security.decoy_analysis import DecoyStateAnalyzer
from snstfqkd.security.bounds import get_S_Z as _get_S_Z


def binary_entropy(x: float) -> float:
    """H(x) = -x log2(x) - (1-x) log2(1-x), with H(0)=H(1)=0."""
    x = float(np.clip(x, 0.0, 1.0))
    if x == 0.0 or x == 1.0:
        return 0.0
    return -x * np.log2(x) - (1.0 - x) * np.log2(1.0 - x)


@dataclass
class KeyRateResult:
    key_rate: float
    epsilon: float
    mu_signal: float
    s1: float
    e1_ph: float
    S_Z: float
    E_Z: float
    distance: float


class SNSKeyRateCalculator:
    """
    Asymptotic SNS-TF-QKD secret key rate per pulse.

    SNS paper Eq. (4):
        R = 2*eps*(1-eps) * mu' * e^{-mu'} * s1 * [1 - H(e1_ph)]
            - S_Z * f * H(E_Z)
    """

    def __init__(self, f: float = 1.1):
        self.f = f

    def key_rate(
        self,
        epsilon: float,
        mu_signal: float,
        s1: float,
        e1_ph: float,
        S_Z: float,
        E_Z: float = 0.0,
    ) -> float:
        """Returns max(0, R) where R is given by SNS Eq. (4)."""
        positive = (
            2.0 * epsilon * (1.0 - epsilon)
            * mu_signal * np.exp(-mu_signal)
            * s1 * (1.0 - binary_entropy(e1_ph))
        )
        correction = S_Z * self.f * binary_entropy(E_Z)
        return max(0.0, positive - correction)

    def key_rate_vs_distance(
        self,
        distances: np.ndarray,
        channel,
        epsilon: float,
        mu_signal: float,
        mu1: float,
        mu2: float,
    ) -> List[KeyRateResult]:
        """
        Sweep key rate over a distance array.

        Uses channel.z_window_yield() (one-arm BS statistics) for S_Z
        rather than channel.intensity_yield() (two-arm interference
        statistics), correctly implementing the Z-window counting rate
        of Wang Eq. (4).
        """
        analyzer = DecoyStateAnalyzer()
        results = []
        for d in distances:
            s0    = channel.vacuum_yield()
            S_mu1 = channel.intensity_yield(mu1, d)       # X-window decoy yield
            S_mu2 = channel.intensity_yield(mu2, d)       # X-window decoy yield
            E_mu1 = channel.qber_x_basis(mu1, d)         # X-window QBER
            E_Z   = channel.qber_z_basis()               # Z-basis QBER (= 0)
            S_Z   = _get_S_Z(channel, mu_signal, d)

            decoy = analyzer.three_intensity_decoy(mu1, mu2, s0, S_mu1, S_mu2, E_mu1)
            rate  = 0.0
            if decoy.valid:
                rate = self.key_rate(epsilon, mu_signal, decoy.s1, decoy.e1_ph, S_Z, E_Z)

            results.append(KeyRateResult(
                key_rate=rate, epsilon=epsilon, mu_signal=mu_signal,
                s1=decoy.s1, e1_ph=decoy.e1_ph,
                S_Z=S_Z, E_Z=E_Z, distance=float(d),
            ))
        return results

    def optimized_key_rate(
        self,
        distance: float,
        channel,
        epsilon: float,
        mu_signal: float,
        mu1: float,
        mu2: float,
    ) -> KeyRateResult:
        return self.key_rate_vs_distance(
            np.array([distance]), channel, epsilon, mu_signal, mu1, mu2
        )[0]
