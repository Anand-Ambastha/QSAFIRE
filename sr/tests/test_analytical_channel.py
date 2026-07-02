"""
test_analytical_channel.py

Verifies AnalyticalChannel reproduces the exact numerical outputs of
the original sns_tfqkd.ChannelModel (golden values captured before
the refactor), and cross-checks basic physical sanity properties.
"""
import numpy as np
import pytest

from snstfqkd.channel.analytical_channel import AnalyticalChannel
from snstfqkd.channel.base import ChannelParameters


# Golden values captured from the ORIGINAL sns_tfqkd.py (pre-refactor),
# params = ChannelParameters(misalignment_error=0.15), distance=300, mu=0.05.
GOLDEN = {
    "vacuum_yield": 1.99999999998e-11,
    "intensity_yield": 7.999600411561013e-05,
    "qber_x_basis": 0.1499965175093047,
    "z_window_yield": 3.999882001742885e-05,
    "gain": 7.999600411561013e-05,
}


@pytest.fixture
def channel():
    return AnalyticalChannel(ChannelParameters(misalignment_error=0.15))


def test_golden_vacuum_yield(channel):
    assert channel.vacuum_yield() == pytest.approx(GOLDEN["vacuum_yield"], rel=1e-12)


def test_golden_intensity_yield(channel):
    assert channel.intensity_yield(0.05, 300) == pytest.approx(GOLDEN["intensity_yield"], rel=1e-12)


def test_golden_qber_x_basis(channel):
    assert channel.qber_x_basis(0.05, 300) == pytest.approx(GOLDEN["qber_x_basis"], rel=1e-12)


def test_golden_z_window_yield(channel):
    assert channel.z_window_yield(0.05, 300) == pytest.approx(GOLDEN["z_window_yield"], rel=1e-12)


def test_golden_gain(channel):
    assert channel.gain(0.05, 300) == pytest.approx(GOLDEN["gain"], rel=1e-12)


def test_transmittance_decreases_with_distance(channel):
    ts = [channel.total_transmittance(d) for d in [0, 100, 300, 600]]
    assert all(ts[i] >= ts[i + 1] for i in range(len(ts) - 1))


def test_intensity_yield_bounds(channel):
    for d in [0, 100, 500, 1000]:
        y = channel.intensity_yield(0.05, d)
        assert 0.0 <= y <= 1.0


def test_qber_z_basis_is_zero(channel):
    assert channel.qber_z_basis() == 0.0


def test_zero_misalignment_reduces_qber():
    ch0 = AnalyticalChannel(ChannelParameters(misalignment_error=0.0))
    ch1 = AnalyticalChannel(ChannelParameters(misalignment_error=0.3))
    assert ch0.qber_x_basis(0.05, 300) <= ch1.qber_x_basis(0.05, 300)
