# -*- coding: utf-8 -*-
"""
snstfqkd.channel.analytical_channel
=====================================

Analytical (closed-form) SNS-TF-QKD channel model.

This is a verbatim move of ``ChannelModel`` from the original
``sns_tfqkd.py``, renamed to ``AnalyticalChannel`` and made to inherit
from ``BaseChannel``. Every equation, every constant, and every line
of arithmetic is unchanged from the original file — only the class
name, the import of shared dataclasses from ``channel.base``, and the
docstrings noting the equation provenance were touched.

Reference: Wang, X.-B., Yu, Z.-W., & Hu, X.-L. (2018). Sending or not
sending: Twin-field quantum key distribution with large misalignment
error. arXiv:1805.09222v9, Eqs. (2), (4), (44), (45).
"""

import numpy as np

from snstfqkd.channel.base import BaseChannel, ChannelParameters, logger


class AnalyticalChannel(BaseChannel):
    """
    Analytical TF-QKD / SNS-TF-QKD channel model at Charlie's node.

    Geometry
    --------
    Alice and Bob each send coherent states of per-mode intensity mu to
    an untrusted relay Charlie located at the midpoint. Each arm has
    transmittance::

        eta_arm(L) = eta_d * 10^(-alpha * L / 20)

    where the factor 20 (not 10) encodes the half-distance geometry.

    Yield model (X-windows) — Level-3 protocol-correct statistics
    --------------------------------------------------------------
    Wang Step 2, Eq. (1) with lambda -> 0 (infinitely small phase slice).
    Only two phase relationships survive, each with probability 1/2:

    Constructive case  cos(delta_A - delta_B) = +1::

        n_cL = 2 eta_arm mu (1 - e_a) + d   [correct detector port]
        n_cR = 2 eta_arm mu   e_a     + d   [error detector port]

    Destructive case   cos(delta_A - delta_B) = -1::

        n_dL = 2 eta_arm mu   e_a     + d   [error detector port]
        n_dR = 2 eta_arm mu (1 - e_a) + d   [correct detector port]

    Single-click yield (averaged over both cases)::

        S_mu = 1/2 [P_cL(1-P_cR) + P_cR(1-P_cL)]
             + 1/2 [P_dL(1-P_dR) + P_dR(1-P_dL)]       [Eq. 44 input]

    Yield model (Z-windows)
    ------------------------
    In a Z-window, exactly one party sends intensity mu' and the other
    sends vacuum. Each BS output receives mean photon number
    ``n_each = eta_arm * mu' / 2 + d``, with no interference (e_a does
    not enter). See ``z_window_yield()``.

    X-basis QBER (Wang Step 5)
    ----------------------------
    Constructive case: error if R fires alone:  err_c = P_cR(1-P_cL)
    Destructive  case: error if L fires alone:  err_d = P_dL(1-P_dR)

        E_mu^X = [1/2*err_c + 1/2*err_d] / S_mu            [Eq. 45 input]

    Z-basis QBER
    ------------
    E_Z = 0 (see docstring of qber_z_basis for justification).
    """

    def __init__(self, params: ChannelParameters) -> None:
        self.params = params
        logger.info("AnalyticalChannel initialised: %s", params)

    # ------------------------------------------------------------------
    # Transmittance
    # ------------------------------------------------------------------

    def channel_transmittance(self, distance: float) -> float:
        """
        Fiber-only transmittance for one arm of length L/2.

        Formula::

            eta_fiber(L) = 10^(-alpha * L / 20)

        The denominator 20 (instead of the conventional 10) reflects the
        TF-QKD half-distance geometry: the arm length is L/2, so the
        exponent is alpha*(L/2)/10 = alpha*L/20.
        """
        return float(np.clip(
            10.0 ** (-self.params.fiber_loss * distance / 20.0),
            0.0, 1.0,
        ))

    def total_transmittance(self, distance: float) -> float:
        """
        Total single-arm transmittance including detector efficiency.

        eta_arm(L) = eta_d * 10^(-alpha * L / 20)
        """
        return self.params.detector_efficiency * self.channel_transmittance(distance)

    # ------------------------------------------------------------------
    # Yield
    # ------------------------------------------------------------------

    def vacuum_yield(self) -> float:
        """
        Vacuum counting rate at Charlie (mu = 0).

        Y_0 = 2 d (1 - d) ~ 2d for d << 1.
        """
        d = self.params.dark_count_rate
        return 2.0 * d * (1.0 - d)

    def intensity_yield(self, intensity: float, distance: float) -> float:
        """
        X-window single-click counting rate S_mu at Charlie's BS.

        Level-3 implementation — protocol-correct X-window statistics
        (Wang 2018, Step 2, Eq. 1, lambda -> 0).

        Returns
        -------
        float
            S_mu (X-window single-click rate) in [0, 1].
            Used as Eq. (44) inputs S_mu1, S_mu2.
            Do NOT use for S_Z in Eq. (4) — use z_window_yield() instead.
        """
        eta_arm = self.total_transmittance(distance)
        d = self.params.dark_count_rate
        ea = self.params.misalignment_error

        # Constructive case: correct port L, wrong port R
        n_cL = 2.0 * eta_arm * intensity * (1.0 - ea) + d
        n_cR = 2.0 * eta_arm * intensity * ea          + d
        # Destructive case: correct port R, wrong port L
        n_dL = 2.0 * eta_arm * intensity * ea          + d
        n_dR = 2.0 * eta_arm * intensity * (1.0 - ea) + d

        P_cL = 1.0 - np.exp(-n_cL)
        P_cR = 1.0 - np.exp(-n_cR)
        P_dL = 1.0 - np.exp(-n_dL)
        P_dR = 1.0 - np.exp(-n_dR)

        S_c = P_cL * (1.0 - P_cR) + P_cR * (1.0 - P_cL)
        S_d = P_dL * (1.0 - P_dR) + P_dR * (1.0 - P_dL)

        yield_ = 0.5 * S_c + 0.5 * S_d
        return float(np.clip(yield_, 0.0, 1.0))

    # ------------------------------------------------------------------
    # Z-window yield  (Wang 2018 Eq. 4, Step 3 definition)
    # ------------------------------------------------------------------

    def z_window_yield(self, intensity: float, distance: float) -> float:
        """
        Z-window single-click counting rate S_Z at Charlie.

        In a Z-window exactly one party sends intensity mu' = intensity;
        the other sends vacuum. Each BS output receives mean photon
        number ``n_each = eta_arm * mu' / 2 + d``. No interference, so
        e_a does not appear; E_Z = 0 regardless of e_a.

        Returns
        -------
        float
            S_Z (Z-window single-click rate) in [0, 1].
        """
        eta_arm = self.total_transmittance(distance)
        d = self.params.dark_count_rate
        # Each BS output receives eta_arm * mu'/2 signal photons + dark counts
        n_each = eta_arm * intensity / 2.0 + d
        P_click = 1.0 - np.exp(-n_each)
        # Single click: one detector fires, the other does not
        S_Z = 2.0 * P_click * (1.0 - P_click)
        return float(np.clip(S_Z, 0.0, 1.0))

    def gain(self, intensity: float, distance: float) -> float:
        """
        Counting rate at Charlie (equals intensity_yield for TF-QKD X-windows).
        """
        return self.intensity_yield(intensity, distance)

    # ------------------------------------------------------------------
    # Error models
    # ------------------------------------------------------------------

    def qber_x_basis(self, intensity: float, distance: float) -> float:
        """
        X-basis quantum bit error rate E_mu^X.

        Level-3 implementation — protocol-correct error counting
        (Wang 2018, Step 5, Eq. 2).

        Returns
        -------
        float
            E_mu^X in [0, 0.5].  Used as Eq. (45) input E_mu1^X.
        """
        eta_arm = self.total_transmittance(distance)
        d = self.params.dark_count_rate
        ea = self.params.misalignment_error

        # Constructive case
        n_cL = 2.0 * eta_arm * intensity * (1.0 - ea) + d
        n_cR = 2.0 * eta_arm * intensity * ea          + d
        # Destructive case
        n_dL = 2.0 * eta_arm * intensity * ea          + d
        n_dR = 2.0 * eta_arm * intensity * (1.0 - ea) + d

        P_cL = 1.0 - np.exp(-n_cL)
        P_cR = 1.0 - np.exp(-n_cR)
        P_dL = 1.0 - np.exp(-n_dL)
        P_dR = 1.0 - np.exp(-n_dR)

        S_c = P_cL * (1.0 - P_cR) + P_cR * (1.0 - P_cL)
        S_d = P_dL * (1.0 - P_dR) + P_dR * (1.0 - P_dL)
        S_mu = 0.5 * S_c + 0.5 * S_d

        if S_mu < 1.0e-30:
            return 0.5

        err_c = P_cR * (1.0 - P_cL)   # R fires alone in constructive = error
        err_d = P_dL * (1.0 - P_dR)   # L fires alone in destructive  = error
        S_err = 0.5 * err_c + 0.5 * err_d

        return float(np.clip(S_err / S_mu, 0.0, 0.5))

    def qber_z_basis(self) -> float:
        """
        Z-basis quantum bit error rate E_Z.

        Returns 0.0 for SNS-TF-QKD. In a Z-window there is no
        alignment-dependent single-photon interference, so no
        misalignment error contributes (Wang 2018, Section III, Eq. 4:
        S_Z f H(E_Z) -> 0).
        """
        return 0.0


# Backward-compatible alias (old code imported ``ChannelModel`` from
# ``sns_tfqkd``). Kept so any external script importing the old name
# during migration keeps working.
ChannelModel = AnalyticalChannel
