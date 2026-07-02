# -*- coding: utf-8 -*-
"""
snstfqkd.protocol.events
===========================

Beam-splitter click model and effective/error event classification,
extracted from ``SNSProtocolSimulator._bs_click_probs`` /
``_sample_clicks`` in the original ``protocol_simulator.py``. Pure
functions — identical math, only pulled out of the class.
"""

import numpy as np


def bs_click_probs(eta_arm: float, ea: float, d: float, mu_A, mu_B, phi):
    """
    Vectorized BS output mean photon numbers for arrays of Alice/Bob
    intensities mu_A, mu_B and relative phase phi = delta_A - delta_B.

    Two-mode coherent-state BS combination (50:50, lossless BS)::

        n_L_raw = eta_arm * (mu_A + mu_B + 2*sqrt(mu_A*mu_B)*cos(phi)) / 2
        n_R_raw = eta_arm * (mu_A + mu_B - 2*sqrt(mu_A*mu_B)*cos(phi)) / 2

    Misalignment error e_a: a fraction e_a of the photons in each port
    are redirected to the opposite port (symmetric crosstalk model)::

        n_L = n_L_raw*(1-e_a) + n_R_raw*e_a + d
        n_R = n_R_raw*(1-e_a) + n_L_raw*e_a + d

    Returns
    -------
    n_L, n_R : ndarray
        Mean photon numbers (incl. dark counts) at left/right ports.
    """
    cross = 2.0 * np.sqrt(np.clip(mu_A * mu_B, 0.0, None)) * np.cos(phi)
    n_L_raw = eta_arm * (mu_A + mu_B + cross) / 2.0
    n_R_raw = eta_arm * (mu_A + mu_B - cross) / 2.0

    n_L = n_L_raw * (1.0 - ea) + n_R_raw * ea + d
    n_R = n_R_raw * (1.0 - ea) + n_L_raw * ea + d
    return n_L, n_R


def sample_clicks(rng: np.random.Generator, n_L, n_R):
    """
    Draw Bernoulli click outcomes for left/right detectors given mean
    photon numbers n_L, n_R (independent Poisson-click model,
    P(click) = 1 - exp(-n)).

    Returns
    -------
    click_L, click_R : boolean ndarray
    """
    P_L = 1.0 - np.exp(-n_L)
    P_R = 1.0 - np.exp(-n_R)
    u_L = rng.random(P_L.shape)
    u_R = rng.random(P_R.shape)
    click_L = u_L < P_L
    click_R = u_R < P_R
    return click_L, click_R


def classify_x_window_events(cos_phi, click_L, click_R):
    """
    Classify X-window click outcomes into effective/error/double/no-click
    events, per Wang Step 5 / Eq. (2): the "correct" detector is L when
    cos(delta_A - delta_B) >= 0 and R when < 0; the other detector
    clicking alone is a wrong X-bit.

    Returns
    -------
    effective, error_event, double_click, no_click : boolean ndarray
    """
    effective = click_L ^ click_R
    double_click = click_L & click_R
    no_click = (~click_L) & (~click_R)

    correct_is_L = cos_phi >= 0.0
    error_event = effective & (
        (correct_is_L & click_R & ~click_L) |
        (~correct_is_L & click_L & ~click_R)
    )
    return effective, error_event, double_click, no_click
