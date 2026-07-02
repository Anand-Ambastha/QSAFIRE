# -*- coding: utf-8 -*-
"""
snstfqkd.channel.detector
===========================

Monte-Carlo coherent-state beam-splitter detector engine used by
``MonteCarloChannel``.

This is a verbatim move of the "PATCH 4" engine
(``_derive_mc_seed`` / ``_mc_beamsplitter_xwindow_stats``) from the
original ``sns_tfqkd_mc.py``. No equations were changed — only moved
into its own module so ``channel/mc_channel.py`` stays focused on the
``BaseChannel`` interface.

Physical model (Steps 1-10 of the observable-generation spec)
----------------------------------------------------------------
  1. delta_phi = phi_A - phi_B, with phi_A, phi_B ~ Uniform(0, 2*pi)
     independently (Alice's / Bob's private random phases).
  2. Channel transmission over each L/2 arm: eta_fiber = 10^(-alpha*L/20).
  3. Coherent-state field amplitudes arriving at Charlie:
     E_A = sqrt(eta_fiber * mu) * e^{i*phi_A}, E_B likewise for Bob.
  4. 50:50 beam-splitter: E_L = (E_A+E_B)/sqrt(2), E_R = (E_A-E_B)/sqrt(2).
  5. Detector-arm intensities: I_L = eta_fiber*mu*(1+cos(delta_phi)),
     I_R = eta_fiber*mu*(1-cos(delta_phi)).
  6. Detector click probabilities including efficiency eta_d and dark
     counts p_d: P = 1 - (1-p_d)*exp(-eta_d * I).
  7. Bernoulli sampling of L/R clicks -> NONE / LEFT / RIGHT / DOUBLE.
  8. Phase-slice postselection, Eq. (1): 1 - |cos(delta_phi)| <= lambda.
     Implemented via EXACT direct sampling of delta_phi from the
     accepted region (not generate-then-reject).
  9. Effective event = exactly one detector clicks (LEFT xor RIGHT).
     Double clicks and no-clicks are excluded from the numerator AND
     from S_mu's denominator (the phase-slice-accepted trial count).
 10. S_mu = N_effective / N_accepted_trials;
     E_mu = N_error_clicks / N_effective (Wang Step 5 / Eq. 2 sign
     rule: wrong-port click is an error).

ASSUMPTION on detector efficiency placement: eta_d is applied ONCE, in
Step 6 (the click-probability step), to avoid unphysical
double-counting between Steps 2-5's fiber-only "eta_arm" and Step 6's
click model.

Misalignment e_a: implemented as symmetric port crosstalk on the
detector-arm intensities.

Reproducibility: the RNG seed is DERIVED deterministically from (mu,
distance, misalignment_error, base_seed) so repeated calls with
identical arguments return IDENTICAL results (required for
``optimizer.py``'s differential_evolution to behave sensibly).
"""

from functools import lru_cache
from typing import Tuple

import numpy as np


def _derive_mc_seed(base_seed: int, *floats: float) -> int:
    """Deterministic seed derived from (base_seed, *floats). Same inputs
    -> same seed -> same MC result (see reproducibility note above)."""
    key = (int(base_seed),) + tuple(round(float(x), 12) for x in floats)
    return hash(key) & 0xFFFFFFFF


@lru_cache(maxsize=200_000)
def mc_beamsplitter_xwindow_stats(
    fiber_loss: float,
    detector_efficiency: float,
    dark_count_rate: float,
    misalignment_error: float,
    intensity: float,
    distance: float,
    n_samples: int,
    lambda_slice: float,
    seed: int,
) -> Tuple[float, float]:
    """
    Monte-Carlo coherent-state beam-splitter simulation of the X-window
    observables (S_mu, E_mu) for one (intensity, distance) operating
    point. Pure function of hashable scalar arguments -> cached with
    lru_cache so repeated evaluations (e.g. optimizer sweeps landing on
    the same point) are O(1).

    Returns
    -------
    (S_mu, E_mu) : (float, float)
        S_mu in [0, 1], E_mu in [0, 0.5] (E_mu = 0.5 fallback if zero
        effective events were observed).
    """
    rng = np.random.default_rng(seed)

    # ---- Step 2: fiber-only channel transmission, L/2 per arm.
    eta_fiber = 10.0 ** (-fiber_loss * distance / 20.0)
    eta_d = detector_efficiency
    p_d = dark_count_rate
    ea = misalignment_error

    # ---- Step 8: direct sampling of delta_phi from the phase-slice
    # accepted region 1-|cos(delta_phi)| <= lambda_slice.
    h = np.arccos(1.0 - lambda_slice)
    branch = rng.integers(0, 2, size=n_samples)
    offset = rng.uniform(-h, h, size=n_samples)
    delta_phi = np.where(branch == 0, offset, np.pi + offset)
    cos_dphi = np.cos(delta_phi)

    # ---- Steps 3-5: coherent-state amplitudes -> BS combination ->
    # detector-arm intensities.
    I_L = eta_fiber * intensity * (1.0 + cos_dphi)
    I_R = eta_fiber * intensity * (1.0 - cos_dphi)

    # Misalignment e_a: symmetric port crosstalk (fraction e_a of each
    # port's intensity leaks to the opposite port).
    I_L_mis = I_L * (1.0 - ea) + I_R * ea
    I_R_mis = I_R * (1.0 - ea) + I_L * ea

    # ---- Step 6: detector click probabilities.
    P_L = 1.0 - (1.0 - p_d) * np.exp(-eta_d * I_L_mis)
    P_R = 1.0 - (1.0 - p_d) * np.exp(-eta_d * I_R_mis)

    # ---- Step 7: Monte-Carlo Bernoulli detector outcomes.
    click_L = rng.random(n_samples) < P_L
    click_R = rng.random(n_samples) < P_R

    # ---- Step 9: effective event = exactly one detector clicks.
    effective = click_L ^ click_R

    # ---- Wang Step 5 / Eq. (2) error rule: correct port is L when
    # cos(delta_phi) >= 0, R when cos(delta_phi) < 0.
    correct_is_L = cos_dphi >= 0.0
    error_event = effective & (
        (correct_is_L & click_R) | (~correct_is_L & click_L)
    )

    n_effective = int(np.count_nonzero(effective))
    n_errors = int(np.count_nonzero(error_event))

    # ---- Step 10: observable construction.
    S_mu = n_effective / n_samples if n_samples > 0 else 0.0
    E_mu = (n_errors / n_effective) if n_effective > 0 else 0.5

    return float(np.clip(S_mu, 0.0, 1.0)), float(np.clip(E_mu, 0.0, 0.5))


def phase_slice_acceptance_fraction(lambda_slice: float) -> float:
    """
    Fraction of the full 2*pi phase circle accepted by the Eq. (1)
    phase-slice criterion 1 - |cos(phi)| <= lambda_slice, for the
    two-band region used by ``sample_sliced_phases`` /
    ``mc_beamsplitter_xwindow_stats``.

    Each of the two symmetric bands has half-width
    h = arccos(1 - lambda_slice); total accepted angular measure is
    2*(2h) = 4h out of 2*pi, giving fraction = 2h/pi.

    This exposes a quantity that was already implicit in the sampling
    method (never previously plotted) — no new physics, just the
    geometric acceptance rate implied by the existing phase-slice
    formula.
    """
    h = np.arccos(1.0 - lambda_slice)
    return float(np.clip(2.0 * h / np.pi, 0.0, 1.0))


def mc_beamsplitter_diagnostics(
    fiber_loss: float,
    detector_efficiency: float,
    dark_count_rate: float,
    misalignment_error: float,
    intensity: float,
    distance: float,
    n_samples: int,
    lambda_slice: float,
    seed: int,
):
    """
    Diagnostic variant of ``mc_beamsplitter_xwindow_stats`` that exposes
    additional intermediate quantities (per-port click rates,
    double-click and no-click fractions) that the production function
    does not return. Uses the IDENTICAL physical model (Steps 1-9,
    reproduced verbatim) — this function performs no new physics, it
    only returns more of what was already being computed internally,
    to support the Monte-Carlo diagnostic figures (MC6-MC11).

    Not lru_cache'd (diagnostic use only, called once per plotted
    distance point rather than from hot optimizer loops).

    Returns
    -------
    dict with keys: S_mu, E_mu, P_L_mean, P_R_mean, effective_fraction,
    double_click_fraction, no_click_fraction, phase_slice_acceptance
    """
    rng = np.random.default_rng(seed)

    eta_fiber = 10.0 ** (-fiber_loss * distance / 20.0)
    eta_d = detector_efficiency
    p_d = dark_count_rate
    ea = misalignment_error

    h = np.arccos(1.0 - lambda_slice)
    branch = rng.integers(0, 2, size=n_samples)
    offset = rng.uniform(-h, h, size=n_samples)
    delta_phi = np.where(branch == 0, offset, np.pi + offset)
    cos_dphi = np.cos(delta_phi)

    I_L = eta_fiber * intensity * (1.0 + cos_dphi)
    I_R = eta_fiber * intensity * (1.0 - cos_dphi)
    I_L_mis = I_L * (1.0 - ea) + I_R * ea
    I_R_mis = I_R * (1.0 - ea) + I_L * ea

    P_L = 1.0 - (1.0 - p_d) * np.exp(-eta_d * I_L_mis)
    P_R = 1.0 - (1.0 - p_d) * np.exp(-eta_d * I_R_mis)

    click_L = rng.random(n_samples) < P_L
    click_R = rng.random(n_samples) < P_R

    effective = click_L ^ click_R
    double_click = click_L & click_R
    no_click = (~click_L) & (~click_R)

    correct_is_L = cos_dphi >= 0.0
    error_event = effective & (
        (correct_is_L & click_R) | (~correct_is_L & click_L)
    )

    n_effective = int(np.count_nonzero(effective))
    n_errors = int(np.count_nonzero(error_event))
    S_mu = n_effective / n_samples if n_samples > 0 else 0.0
    E_mu = (n_errors / n_effective) if n_effective > 0 else 0.5

    return {
        "S_mu": float(np.clip(S_mu, 0.0, 1.0)),
        "E_mu": float(np.clip(E_mu, 0.0, 0.5)),
        "P_L_mean": float(np.mean(P_L)),
        "P_R_mean": float(np.mean(P_R)),
        "effective_fraction": float(np.mean(effective)),
        "double_click_fraction": float(np.mean(double_click)),
        "no_click_fraction": float(np.mean(no_click)),
        "phase_slice_acceptance": phase_slice_acceptance_fraction(lambda_slice),
    }
