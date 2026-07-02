# -*- coding: utf-8 -*-
"""
snstfqkd.channel.mc_channel
=============================

Monte-Carlo (beam-splitter detector-event simulation) SNS-TF-QKD
channel model — "Layer 2" of the audit.

This is a verbatim move of ``ChannelModel`` from the original
``sns_tfqkd_mc.py``, renamed to ``MonteCarloChannel`` and made to
inherit from ``BaseChannel``. The MC engine itself
(``_mc_beamsplitter_xwindow_stats`` / ``_derive_mc_seed``) now lives in
``channel/detector.py``; everything else (transmittance, vacuum yield,
z_window_yield, gain, qber_z_basis) is identical to
``AnalyticalChannel`` since those methods were identical in the two
original files.

IMPORTANT — kept deliberately distinct from ``protocol.SNSProtocolSimulator``
--------------------------------------------------------------------------
Per the audit finding: ``MonteCarloChannel`` here implements the SAME
``BaseChannel`` interface as ``AnalyticalChannel`` (drop-in replacement,
fast, cached, per-point sampling). ``SNSProtocolSimulator`` in
``protocol/`` is a different, lower-level abstraction — a literal
pulse-by-pulse protocol simulator with its own API (window selection,
send/no-send decisions, phase generation) — and is NOT merged into
this class. Both are preserved independently, as required.

Reference: Wang, X.-B., Yu, Z.-W., & Hu, X.-L. (2018), arXiv:1805.09222v9.
"""

import numpy as np

from snstfqkd.channel.base import BaseChannel, ChannelParameters, logger
from snstfqkd.channel.detector import _derive_mc_seed, mc_beamsplitter_xwindow_stats


class MonteCarloChannel(BaseChannel):
    """
    Monte-Carlo coherent-state beam-splitter SNS-TF-QKD channel model.

    ``intensity_yield()`` and ``qber_x_basis()`` derive S_mu / E_mu^X
    from explicit Monte-Carlo event counting over simulated coherent
    -state pulse pairs at Charlie (see ``channel.detector`` for the
    full physical model), rather than the closed-form expressions used
    by ``AnalyticalChannel``. All other methods (transmittance,
    vacuum_yield, z_window_yield, gain, qber_z_basis) use the identical
    formulas as ``AnalyticalChannel`` — these were never patched to be
    Monte-Carlo based in the original codebase.

    In the lambda -> 0, N_samples -> infinity limit this Monte-Carlo
    estimator converges to the same closed-form expression used by
    ``AnalyticalChannel``.
    """

    def __init__(
        self,
        params: ChannelParameters,
        mc_n_samples: int = 2_000_000,
        #mc_n_samples: int = 100_000_000, #for publication
        mc_lambda_slice: float = 1.0e-3,
        mc_seed: int = 20260701,
    ) -> None:
        """
        Parameters
        ----------
        params : ChannelParameters
            Physical channel/detector parameters (unchanged).
        mc_n_samples : int, optional
            Number of phase-slice-accepted Monte-Carlo trials used
            internally by intensity_yield()/qber_x_basis() for each
            (mu, distance) evaluation. Default 2,000,000.
        mc_lambda_slice : float, optional
            Finite phase-slice half-width used by the direct sampler
            (Eq. 1). Default 1e-3.
        mc_seed : int, optional
            Base seed for the deterministic per-call MC seed
            derivation. Default 20260701.
        """
        self.params = params
        self._mc_n_samples = int(mc_n_samples)
        self._mc_lambda_slice = float(mc_lambda_slice)
        self._mc_seed = int(mc_seed)
        logger.info("MonteCarloChannel initialised: %s", params)

    # ------------------------------------------------------------------
    # Transmittance  (identical to AnalyticalChannel)
    # ------------------------------------------------------------------

    def channel_transmittance(self, distance: float) -> float:
        """Fiber-only transmittance for one arm of length L/2."""
        return float(np.clip(
            10.0 ** (-self.params.fiber_loss * distance / 20.0),
            0.0, 1.0,
        ))

    def total_transmittance(self, distance: float) -> float:
        """Total single-arm transmittance including detector efficiency."""
        return self.params.detector_efficiency * self.channel_transmittance(distance)

    # ------------------------------------------------------------------
    # Yield
    # ------------------------------------------------------------------

    def vacuum_yield(self) -> float:
        """Vacuum counting rate at Charlie (mu = 0). Y_0 = 2 d (1 - d)."""
        d = self.params.dark_count_rate
        return 2.0 * d * (1.0 - d)

    def intensity_yield(self, intensity: float, distance: float) -> float:
        """
        X-window single-click counting rate S_mu at Charlie's BS,
        derived from Monte-Carlo event counting.

        Returns
        -------
        float
            S_mu (X-window single-click rate) in [0, 1].
            Used as Eq. (44) inputs S_mu1, S_mu2.
            Do NOT use for S_Z in Eq. (4) — use z_window_yield() instead.
        """
        S_mu, _E_mu = self._mc_x_window_observables(intensity, distance)
        return S_mu

    # ------------------------------------------------------------------
    # Z-window yield  (identical to AnalyticalChannel — never patched
    # to be Monte-Carlo based in the original codebase)
    # ------------------------------------------------------------------

    def z_window_yield(self, intensity: float, distance: float) -> float:
        """Z-window single-click counting rate S_Z at Charlie."""
        eta_arm = self.total_transmittance(distance)
        d = self.params.dark_count_rate
        n_each = eta_arm * intensity / 2.0 + d
        P_click = 1.0 - np.exp(-n_each)
        S_Z = 2.0 * P_click * (1.0 - P_click)
        return float(np.clip(S_Z, 0.0, 1.0))

    def gain(self, intensity: float, distance: float) -> float:
        """Counting rate at Charlie (equals intensity_yield for X-windows)."""
        return self.intensity_yield(intensity, distance)

    # ------------------------------------------------------------------
    # Error models
    # ------------------------------------------------------------------

    def qber_x_basis(self, intensity: float, distance: float) -> float:
        """
        X-basis quantum bit error rate E_mu^X, derived from the SAME
        simulated pulse stream as intensity_yield() (shared seed via
        ``_mc_x_window_observables``).

        Returns
        -------
        float
            E_mu^X in [0, 0.5]. Used as Eq. (45) input E_mu1^X.
        """
        _S_mu, E_mu = self._mc_x_window_observables(intensity, distance)
        return E_mu

    # ------------------------------------------------------------------
    # Shared Monte-Carlo entry point for intensity_yield() and
    # qber_x_basis(). Both public methods call this SAME helper so that
    # S_mu and E_mu are computed from one consistent simulated pulse
    # stream (same RNG seed for a given (intensity, distance) point).
    # ------------------------------------------------------------------
    def _mc_x_window_observables(self, intensity: float, distance: float):
        seed = _derive_mc_seed(
            self._mc_seed,
            intensity,
            distance,
            self.params.misalignment_error,
        )
        return mc_beamsplitter_xwindow_stats(
            self.params.fiber_loss,
            self.params.detector_efficiency,
            self.params.dark_count_rate,
            self.params.misalignment_error,
            float(intensity),
            float(distance),
            self._mc_n_samples,
            self._mc_lambda_slice,
            seed,
        )

    def qber_z_basis(self) -> float:
        """
        Z-basis quantum bit error rate E_Z. Returns 0.0 (identical
        justification to AnalyticalChannel.qber_z_basis).
        """
        return 0.0


# Backward-compatible alias (old code imported ``ChannelModel`` from
# ``sns_tfqkd_mc``).
ChannelModel = MonteCarloChannel
