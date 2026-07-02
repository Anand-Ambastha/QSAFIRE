# -*- coding: utf-8 -*-
"""
snstfqkd.protocol.pulse
==========================

Alice/Bob send/no-send decision sampling for Z-windows and the
resulting single-arm photon-number statistic, extracted from
``SNSProtocolSimulator.simulate_z_windows`` in the original
``protocol_simulator.py``. Pure functions — identical math, only
pulled out of the class.
"""

import numpy as np


def sample_z_window_intensities(rng: np.random.Generator, batch: int, epsilon: float, mu_signal: float):
    """
    Alice and Bob each independently decide to send (prob eps) in this
    already-committed signal (Z-)window (Wang Step 1).

    Returns
    -------
    mu_A, mu_B : ndarray
        Per-trial sent intensity (mu_signal or 0.0) for Alice and Bob.
    """
    alice_sends = rng.random(batch) < epsilon
    bob_sends = rng.random(batch) < epsilon
    mu_A = np.where(alice_sends, mu_signal, 0.0)
    mu_B = np.where(bob_sends, mu_signal, 0.0)
    return mu_A, mu_B


def z_window_photon_number(eta_arm: float, d: float, mu_A, mu_B):
    """
    Z-windows have no two-arm interference (Wang: no single-photon
    interference required in Z basis), so the relevant BS statistic is
    the INCOHERENT sum of mean photon numbers split equally between
    output ports, exactly as in AnalyticalChannel.z_window_yield().
    Each party's photon goes to BOTH ports equally (50:50 BS on a
    single-arm input has no phase dependence), so no coherent cross
    term appears.
    """
    return eta_arm * (mu_A + mu_B) / 2.0 + d
