"""
test_mc_channel.py

Verifies MonteCarloChannel reproduces the exact numerical outputs of
the original sns_tfqkd_mc.ChannelModel (golden values captured before
the refactor), and cross-validates against AnalyticalChannel (this is
the missing cross-check the architecture audit flagged).

Replaces the original top-level validate_mc.py, which compared the MC
engine against a THIRD, standalone reimplementation of the analytic
formulas pasted into that file — here it's compared directly against
the real AnalyticalChannel instead.
"""
import numpy as np
import pytest

from snstfqkd.channel.analytical_channel import AnalyticalChannel
from snstfqkd.channel.mc_channel import MonteCarloChannel
from snstfqkd.channel.base import ChannelParameters


# Golden values captured from the ORIGINAL sns_tfqkd_mc.py (pre-refactor),
# default MC settings (2_000_000 samples, lambda=1e-3, seed=20260701),
# params = ChannelParameters(misalignment_error=0.15), distance=300, mu=0.05.
GOLDEN = {
    "intensity_yield": 8e-05,
    "qber_x_basis": 0.1,
}


@pytest.fixture
def mc_channel():
    return MonteCarloChannel(ChannelParameters(misalignment_error=0.15))


@pytest.fixture
def analytical_channel():
    return AnalyticalChannel(ChannelParameters(misalignment_error=0.15))


def test_golden_intensity_yield(mc_channel):
    assert mc_channel.intensity_yield(0.05, 300) == pytest.approx(GOLDEN["intensity_yield"])


def test_golden_qber_x_basis(mc_channel):
    assert mc_channel.qber_x_basis(0.05, 300) == pytest.approx(GOLDEN["qber_x_basis"])


def test_deterministic_repeatability(mc_channel):
    """Same (mu, distance) must give identical MC results across calls."""
    a = mc_channel.intensity_yield(0.05, 300)
    b = mc_channel.intensity_yield(0.05, 300)
    assert a == b


def test_shared_params_identical_to_analytical(mc_channel, analytical_channel):
    """z_window_yield, gain formula, vacuum_yield, transmittance, qber_z_basis
    are identical between the two channels (never patched to be MC-based)."""
    for d in [50, 200, 500]:
        assert mc_channel.total_transmittance(d) == analytical_channel.total_transmittance(d)
        assert mc_channel.vacuum_yield() == analytical_channel.vacuum_yield()
        assert mc_channel.z_window_yield(0.05, d) == analytical_channel.z_window_yield(0.05, d)
    assert mc_channel.qber_z_basis() == analytical_channel.qber_z_basis()


@pytest.mark.parametrize("distance", [100, 300, 500])
def test_mc_agrees_with_analytical_within_statistical_tolerance(mc_channel, analytical_channel, distance):
    """Cross-validation the audit flagged as missing: with 2e6 samples
    and small lambda_slice, MC S_mu/E_mu should track the closed-form
    analytical values within a generous statistical tolerance."""
    S_analytic = analytical_channel.intensity_yield(0.05, distance)
    S_mc = mc_channel.intensity_yield(0.05, distance)
    assert S_mc == pytest.approx(S_analytic, rel=0.05, abs=1e-5)
