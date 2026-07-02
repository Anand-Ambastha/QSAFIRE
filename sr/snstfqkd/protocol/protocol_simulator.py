# -*- coding: utf-8 -*-
"""
snstfqkd.protocol.protocol_simulator
=======================================

Protocol-level SNS-TF-QKD Monte Carlo simulator — "Layer 3" of the
architecture audit.

PRESERVED AS AN INDEPENDENT ABSTRACTION (per explicit requirement):
this is NOT the same thing as ``channel.MonteCarloChannel``. It
replaces the analytical / closed-form X-window and Z-window
observables with observables that EMERGE from explicit Monte Carlo
simulation of the protocol's random variables (window choices,
send/no-send decisions, intensity choices, phases, and Charlie's
beam-splitter interference and detector clicks).

This module only orchestrates; the pure numerical pieces now live in
``protocol.phase_slice``, ``protocol.events``, and ``protocol.pulse``
(extracted verbatim, no equations changed) and the result containers
live in ``protocol.statistics``.

Eq. (44), Eq. (45), and Eq. (4) are NOT modified anywhere in this
package. This module only produces the five observed statistics
S_mu0, S_mu1, S_mu2, E_X_mu1, S_Z, which are then handed unmodified to
``security.DecoyStateAnalyzer.three_intensity_decoy()`` [Eq. 44, 45]
and ``security.SNSKeyRateCalculator.key_rate()`` [Eq. 4].
"""

import numpy as np
from typing import Optional

from snstfqkd.channel.analytical_channel import AnalyticalChannel
from snstfqkd.channel.base import ChannelParameters
from snstfqkd.protocol.statistics import XWindowStats, ZWindowStats, ProtocolStatistics
from snstfqkd.protocol.phase_slice import sample_sliced_phases
from snstfqkd.protocol.events import bs_click_probs, sample_clicks, classify_x_window_events
from snstfqkd.protocol.pulse import sample_z_window_intensities, z_window_photon_number


class SNSProtocolSimulator:
    """
    Protocol-level Monte Carlo simulator for SNS-TF-QKD (Wang et al. 2018).

    For each simulated pulse pair, Alice and Bob independently:
        1. choose signal window or decoy window,
        2. if signal window: send (prob eps) or not-send (prob 1-eps),
        3. if decoy window: choose an intensity from {0, mu1, mu2},
        4. draw a private random phase delta ~ Uniform(0, 2*pi)
           (phase is irrelevant for vacuum / not-sending pulses).

    Misalignment error e_a is applied as a coherent mixing between the
    two output ports. This reduces, in the lambda -> 0 limit, to the
    closed-form Level-3 formulas implemented in
    ``AnalyticalChannel.intensity_yield()`` / ``qber_x_basis()`` — the
    Monte Carlo result is validated against those in the test suite.

    Public methods
    --------------
    simulate_x_windows(intensity, n_pulses) -> XWindowStats
    simulate_z_windows(n_pulses)            -> ZWindowStats
    collect_statistics()                    -> ProtocolStatistics
    run()                                   -> ProtocolStatistics
    """

    def __init__(
        self,
        channel,
        distance: float,
        epsilon: float = 0.05,
        mu_signal: float = 0.05,
        mu1: float = 0.01,
        mu2: float = 0.15,
        lambda_slice: float = 1e-4,
        n_pulses_x: int = 50_000_000,
        n_pulses_z: int = 50_000_000,
        seed: int = 42,
    ):
        """
        Parameters
        ----------
        channel : BaseChannel
            Only its physical PARAMETERS (fiber loss, detector
            efficiency, dark count rate, misalignment error) are reused
            here. Its intensity_yield/qber_x_basis/z_window_yield
            methods are deliberately NOT called by this simulator.
        distance : float
            Total Alice-Bob distance L [km].
        epsilon : float
            Signal-window sending probability eps (Wang Step 1).
        mu_signal : float
            Signal-pulse intensity mu' (Wang Step 1).
        mu1, mu2 : float
            Decoy intensities (mu1 < mu2), with vacuum (mu0=0) implicit.
        lambda_slice : float
            Phase-slice half-width lambda in Eq. (1).
        n_pulses_x : int
            Number of (Alice, Bob) decoy-window trial pairs simulated
            per intensity.
        n_pulses_z : int
            Number of (Alice, Bob) signal-window trial pairs simulated.
        seed : int
            RNG seed for reproducibility.
        """
        self.channel = channel
        self.params = channel.params
        self.distance = float(distance)
        self.epsilon = float(epsilon)
        self.mu_signal = float(mu_signal)
        self.mu1 = float(mu1)
        self.mu2 = float(mu2)
        self.lambda_slice = float(lambda_slice)
        self.n_pulses_x = int(n_pulses_x)
        self.n_pulses_z = int(n_pulses_z)
        self.seed = int(seed)
        self._rng = np.random.default_rng(seed)

        # Physical constants pulled from the channel parameters.
        self.eta_arm = channel.total_transmittance(self.distance)
        self.d = self.params.dark_count_rate
        self.ea = self.params.misalignment_error

    # -----------------------------------------------------------------
    # X-window simulation (Wang Step 1 decoy branch + Step 2 + Eq. 1)
    # -----------------------------------------------------------------
    def simulate_x_windows(self, intensity: float, n_pulses: Optional[int] = None) -> XWindowStats:
        """
        Simulate explicit X-window trials at a single decoy intensity.

        Processed in fixed-size chunks (2,000,000 pulses per chunk) to
        bound peak memory regardless of n_pulses; the resulting counts
        are accumulated exactly as if processed in one single batch.
        """
        if n_pulses is None:
            n_pulses = self.n_pulses_x

        CHUNK = 2_000_000
        n_acc_total = 0
        n_effective = 0
        n_errors = 0
        n_double = 0
        n_no = 0

        remaining = n_pulses
        while remaining > 0:
            batch = min(CHUNK, remaining)
            remaining -= batch

            phi_acc = sample_sliced_phases(self._rng, batch, self.lambda_slice)
            n_acc = phi_acc.shape[0]
            if intensity <= 0.0:
                mu_A = np.zeros(n_acc)
                mu_B = np.zeros(n_acc)
            else:
                mu_A = np.full(n_acc, intensity)
                mu_B = np.full(n_acc, intensity)

            n_L, n_R = bs_click_probs(self.eta_arm, self.ea, self.d, mu_A, mu_B, phi_acc)
            click_L, click_R = sample_clicks(self._rng, n_L, n_R)

            cos_phi = np.cos(phi_acc)
            effective, error_event, double_click, no_click = classify_x_window_events(
                cos_phi, click_L, click_R
            )

            n_acc_total += n_acc
            n_effective += int(np.count_nonzero(effective))
            n_errors += int(np.count_nonzero(error_event))
            n_double += int(np.count_nonzero(double_click))
            n_no += int(np.count_nonzero(no_click))

        if n_acc_total == 0:
            return XWindowStats(
                intensity=intensity, n_trials=n_pulses, n_phase_slice=0,
                n_effective=0, n_errors=0, n_double_click=0, n_no_click=0,
                S_mu=0.0, E_mu=0.5,
            )

        S_mu = n_effective / n_acc_total
        E_mu = n_errors / n_effective if n_effective > 0 else 0.5

        return XWindowStats(
            intensity=intensity,
            n_trials=n_pulses,
            n_phase_slice=n_acc_total,
            n_effective=n_effective,
            n_errors=n_errors,
            n_double_click=n_double,
            n_no_click=n_no,
            S_mu=float(S_mu),
            E_mu=float(np.clip(E_mu, 0.0, 0.5)),
        )

    # -----------------------------------------------------------------
    # Z-window simulation (Wang Step 1 signal branch + Step 2)
    # -----------------------------------------------------------------
    def simulate_z_windows(self, n_pulses: Optional[int] = None) -> ZWindowStats:
        """
        Simulate explicit Z-window trials: both parties have already
        committed to a signal window; each independently sends
        mu_signal (prob eps) or nothing.

        Note: this generates the TRUE eps-weighted mixture over all
        four send/no-send combinations, which differs from
        ``AnalyticalChannel.z_window_yield()`` (which implicitly
        assumes exactly one party always sends) — a genuine,
        protocol-derived discrepancy, not Monte Carlo noise.
        """
        if n_pulses is None:
            n_pulses = self.n_pulses_z

        CHUNK = 2_000_000
        n_effective = 0
        remaining = n_pulses
        while remaining > 0:
            batch = min(CHUNK, remaining)
            remaining -= batch

            mu_A, mu_B = sample_z_window_intensities(self._rng, batch, self.epsilon, self.mu_signal)
            n_each = z_window_photon_number(self.eta_arm, self.d, mu_A, mu_B)

            P_click = 1.0 - np.exp(-n_each)
            u_L = self._rng.random(batch)
            u_R = self._rng.random(batch)
            click_L = u_L < P_click
            click_R = u_R < P_click

            effective = click_L ^ click_R
            n_effective += int(np.count_nonzero(effective))

        n_z = n_pulses  # all trials are Z-windows by construction
        S_Z = n_effective / n_z if n_z > 0 else 0.0

        return ZWindowStats(
            mu_signal=self.mu_signal,
            n_trials=n_pulses,
            n_z_windows=n_z,
            n_effective=n_effective,
            S_Z=float(np.clip(S_Z, 0.0, 1.0)),
        )

    # -----------------------------------------------------------------
    # Aggregation
    # -----------------------------------------------------------------
    def collect_statistics(
        self,
        x_stats_mu0: XWindowStats,
        x_stats_mu1: XWindowStats,
        x_stats_mu2: XWindowStats,
        z_stats: ZWindowStats,
    ) -> ProtocolStatistics:
        """
        Package raw simulate_x_windows()/simulate_z_windows() outputs
        into the ProtocolStatistics container, in exactly the field
        layout consumed by DecoyStateAnalyzer.three_intensity_decoy()
        and SNSKeyRateCalculator.key_rate(). No estimation or smoothing
        is performed here.
        """
        return ProtocolStatistics(
            distance=self.distance,
            mu1=self.mu1,
            mu2=self.mu2,
            mu_signal=self.mu_signal,
            epsilon=self.epsilon,
            lambda_slice=self.lambda_slice,
            seed=self.seed,
            S_mu0=x_stats_mu0.S_mu,
            S_mu1=x_stats_mu1.S_mu,
            S_mu2=x_stats_mu2.S_mu,
            E_X_mu1=x_stats_mu1.E_mu,
            S_Z=z_stats.S_Z,
            E_Z=0.0,   # Wang Sec. III: no misalignment error in Z-basis
            x_window_details={
                0.0: x_stats_mu0,
                self.mu1: x_stats_mu1,
                self.mu2: x_stats_mu2,
            },
            z_window_details=z_stats,
        )

    def run(self) -> ProtocolStatistics:
        """
        Run the full protocol-level Monte Carlo simulation at the
        configured (distance, epsilon, mu_signal, mu1, mu2,
        lambda_slice) and return the observed ProtocolStatistics.
        """
        x0 = self.simulate_x_windows(0.0)
        x1 = self.simulate_x_windows(self.mu1)
        x2 = self.simulate_x_windows(self.mu2)
        z = self.simulate_z_windows()
        return self.collect_statistics(x0, x1, x2, z)


# ---------------------------------------------------------------------
# Convenience adapter: feed ProtocolStatistics into the EXISTING,
# UNMODIFIED DecoyStateAnalyzer / SNSKeyRateCalculator.
# ---------------------------------------------------------------------
def evaluate_key_rate_from_simulation(
    stats: ProtocolStatistics,
    analyzer,      # security.decoy_analysis.DecoyStateAnalyzer instance
    calculator,    # security.key_rate.SNSKeyRateCalculator instance
):
    """
    Glue function: takes Monte-Carlo-observed ProtocolStatistics and
    feeds them, completely unmodified, into the existing Eq. (44)/(45)
    decoy analysis and Eq. (4) key-rate calculator. Performs NO physics
    itself; it only routes the five observed numbers plus the protocol
    parameters (epsilon, mu_signal) to the existing, unchanged APIs.

    Returns
    -------
    (decoy_result, key_rate) : (DecoyResult, float)
    """
    decoy_result = analyzer.three_intensity_decoy(
        stats.mu1, stats.mu2,
        stats.S_mu0, stats.S_mu1, stats.S_mu2, stats.E_X_mu1,
    )
    if not decoy_result.valid:
        return decoy_result, 0.0

    rate = calculator.key_rate(
        stats.epsilon, stats.mu_signal,
        decoy_result.s1, decoy_result.e1_ph,
        stats.S_Z, stats.E_Z,
    )
    return decoy_result, rate
