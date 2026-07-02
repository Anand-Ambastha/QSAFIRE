"""
test_optimizer.py

Verifies SNSOptimizer.optimize_all() reproduces the exact numerical
outputs of the original optimizer.py (golden values captured before
the refactor).
"""
import pytest

from snstfqkd.channel.analytical_channel import AnalyticalChannel
from snstfqkd.channel.base import ChannelParameters
from snstfqkd.optimization.optimizer import SNSOptimizer


# Golden values captured from the ORIGINAL optimizer.py (pre-refactor),
# params = ChannelParameters(misalignment_error=0.15), distance=300,
# de_maxiter=60, seed=42.
GOLDEN = dict(
    key_rate=4.7318264457779755e-05,
    epsilon=0.49988551779813795,
    mu_signal=0.5,
    mu1=0.0001,
    mu2=0.07160829520789341,
    s1=0.0007999880071946559,
    e1_ph=0.1500322427598695,
)


def test_optimize_all_golden_values():
    ch = AnalyticalChannel(ChannelParameters(misalignment_error=0.15))
    opt = SNSOptimizer(ch, de_maxiter=60)
    r = opt.optimize_all(300, de_maxiter=60, seed=42)

    assert r.key_rate == pytest.approx(GOLDEN["key_rate"], rel=1e-9)
    assert r.epsilon == pytest.approx(GOLDEN["epsilon"], rel=1e-9)
    assert r.mu_signal == pytest.approx(GOLDEN["mu_signal"], rel=1e-9)
    assert r.mu1 == pytest.approx(GOLDEN["mu1"], rel=1e-9)
    assert r.mu2 == pytest.approx(GOLDEN["mu2"], rel=1e-9)
    assert r.s1 == pytest.approx(GOLDEN["s1"], rel=1e-9)
    assert r.e1_ph == pytest.approx(GOLDEN["e1_ph"], rel=1e-9)


def test_optimize_all_key_rate_non_negative():
    ch = AnalyticalChannel(ChannelParameters(misalignment_error=0.15))
    opt = SNSOptimizer(ch, de_maxiter=30)
    for d in [50, 300, 700]:
        r = opt.optimize_all(d, de_maxiter=30, seed=1)
        assert r.key_rate >= 0.0


def test_optimize_all_respects_mu1_lt_mu2():
    ch = AnalyticalChannel(ChannelParameters(misalignment_error=0.1))
    opt = SNSOptimizer(ch, de_maxiter=30)
    r = opt.optimize_all(300, de_maxiter=30, seed=7)
    assert r.mu1 < r.mu2
