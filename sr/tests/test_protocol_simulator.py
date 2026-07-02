"""
test_protocol_simulator.py

Verifies SNSProtocolSimulator (Layer 3, kept independent of both
channel implementations) converges to AnalyticalChannel's closed-form
values, and that the phase-slice direct sampler matches naive
rejection sampling. Adapted from the original protocol_test.py.
"""
import numpy as np
import pytest

from snstfqkd.channel.analytical_channel import AnalyticalChannel
from snstfqkd.channel.base import ChannelParameters
from snstfqkd.protocol.protocol_simulator import SNSProtocolSimulator, evaluate_key_rate_from_simulation
from snstfqkd.protocol.phase_slice import sample_sliced_phases
from snstfqkd.security.decoy_analysis import DecoyStateAnalyzer
from snstfqkd.security.key_rate import SNSKeyRateCalculator


def test_phase_slice_sampler_matches_rejection():
    rng_direct = np.random.default_rng(1)
    lambda_slice = 0.05

    rng = np.random.default_rng(7)
    N_raw = 2_000_000
    phi_raw = rng.uniform(0, 2 * np.pi, N_raw)
    mask = (1 - np.abs(np.cos(phi_raw))) <= lambda_slice
    phi_rejected = phi_raw[mask]

    phi_direct = sample_sliced_phases(rng_direct, len(phi_rejected), lambda_slice)

    c_rej = np.abs(np.cos(phi_rejected))
    c_dir = np.abs(np.cos(phi_direct))
    rel_diff = abs(c_rej.mean() - c_dir.mean()) / c_rej.mean()
    assert rel_diff < 0.02


def test_protocol_s_mu1_converges_to_analytical():
    params = ChannelParameters(misalignment_error=0.15)
    channel = AnalyticalChannel(params)
    sim = SNSProtocolSimulator(
        channel, distance=300.0, mu1=0.01, mu2=0.15, mu_signal=0.05,
        epsilon=0.05, lambda_slice=1e-3, n_pulses_x=3_000_000, n_pulses_z=1_000_000, seed=1,
    )
    stats = sim.simulate_x_windows(0.01)
    analytic = channel.intensity_yield(0.01, 300.0)
    assert stats.S_mu == pytest.approx(analytic, rel=0.1, abs=1e-5)


def test_evaluate_key_rate_from_simulation_runs():
    params = ChannelParameters(misalignment_error=0.1)
    channel = AnalyticalChannel(params)
    sim = SNSProtocolSimulator(
        channel, distance=200.0, mu1=0.01, mu2=0.15, mu_signal=0.05,
        epsilon=0.05, n_pulses_x=500_000, n_pulses_z=500_000, seed=2,
    )
    stats = sim.run()
    analyzer = DecoyStateAnalyzer()
    calc = SNSKeyRateCalculator()
    decoy, rate = evaluate_key_rate_from_simulation(stats, analyzer, calc)
    assert rate >= 0.0
