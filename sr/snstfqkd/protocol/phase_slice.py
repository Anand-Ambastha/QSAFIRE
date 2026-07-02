# -*- coding: utf-8 -*-
"""
snstfqkd.protocol.phase_slice
================================

Direct (rejection-free) sampling of the Eq. (1) phase-slice region,
extracted from ``SNSProtocolSimulator._sample_sliced_phases`` in the
original ``protocol_simulator.py``. Pure function of an explicit RNG —
identical math, only pulled out of the class so it can be unit tested
and reused independently.
"""

import numpy as np


def sample_sliced_phases(rng: np.random.Generator, n_pulses: int, lambda_slice: float) -> np.ndarray:
    """
    Draw n_pulses phase differences phi = delta_A - delta_B directly
    from the distribution induced by:
        (a) delta_A, delta_B ~ Uniform(0, 2*pi) independently, and
        (b) Eq. (1) post-selection: 1 - |cos(phi)| <= lambda_slice.

    Rather than drawing delta_A, delta_B uniformly and REJECTING the
    ~(1 - P_accept) fraction that fails Eq. (1) -- which wastes
    O(1/lambda_slice) trials for small lambda_slice -- we sample
    directly from the (exact, closed-form) conditional distribution
    of phi given that it lies in the accepted region. This is
    ordinary importance/inverse-transform sampling restricted to the
    support defined by Eq. (1); it changes neither the protocol nor
    the underlying phase statistics, only the sampling EFFICIENCY.

    The accepted region |cos(phi)| >= 1 - lambda_slice consists of
    two symmetric bands:
        band 0 (constructive, cos(phi) > 0): phi in [-h, +h]
        band 1 (destructive,  cos(phi) < 0): phi in [pi-h, pi+h]
    where h = arccos(1 - lambda_slice). For phi uniform on the full
    circle, EACH band carries exactly the same probability mass by
    symmetry, and WITHIN each band phi is uniform. We therefore draw
    the band uniformly at random (prob 1/2 each, matching the true
    conditional probabilities) and then phi uniformly within the
    chosen band's half-width h. This reproduces exactly the same
    conditional distribution that rejection sampling would yield
    (verified in the test suite against direct rejection sampling).
    """
    h = np.arccos(1.0 - lambda_slice)
    branch = rng.integers(0, 2, size=n_pulses)
    offset = rng.uniform(-h, h, size=n_pulses)
    phi = np.where(branch == 0, offset, np.pi + offset)
    return phi
