# tests/test_phase23.py
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pytest

from snstfqkd.security.decoy_analysis import DecoyStateAnalyzer
from snstfqkd.security.key_rate import binary_entropy, SNSKeyRateCalculator


class FakeChannel:
    def __init__(self, eta=1e-4, dark=1e-11, ea=0.0):
        self.eta = eta; self.dark = dark; self.ea = ea
    def vacuum_yield(self):
        return 2.0 * self.dark
    def intensity_yield(self, mu, distance):
        return 1 - (1 - self.vacuum_yield()) * np.exp(-self.eta * mu)
    def qber_x_basis(self, mu, distance):
        return self.ea + (1 - self.ea) * (1 - np.exp(-self.eta * mu)) / 2.0
    def qber_z_basis(self):
        return 0.0


class TestPhotonProbabilities:
    def setup_method(self): self.da = DecoyStateAnalyzer()
    def test_p0_at_zero(self): assert self.da.p0(0.0) == pytest.approx(1.0)
    def test_p0_positive(self):
        for mu in [0.01, 0.1, 0.5]: assert self.da.p0(mu) > 0
    def test_p1_peak(self):
        mus = np.linspace(0.01, 2.0, 200)
        assert abs(mus[np.argmax([self.da.p1(m) for m in mus])] - 0.5) < 0.05
    def test_p2_at_zero(self): assert self.da.p2(0.0) == pytest.approx(0.0)
    def test_pk_sums_to_one(self):
        for mu in [0.1, 0.3, 0.5]:
            assert sum(self.da.pk(mu, k) for k in range(31)) == pytest.approx(1.0, abs=1e-6)
    def test_p0_p1_p2_consistent_with_pk(self):
        for mu in [0.1, 0.25, 0.5]:
            assert self.da.p0(mu) == pytest.approx(self.da.pk(mu, 0), rel=1e-10)
            assert self.da.p1(mu) == pytest.approx(self.da.pk(mu, 1), rel=1e-10)
            assert self.da.p2(mu) == pytest.approx(self.da.pk(mu, 2), rel=1e-10)


class TestDecoyBounds:
    def setup_method(self):
        self.da = DecoyStateAnalyzer()
        self.ch = FakeChannel(eta=1e-3)
    def _decoy(self, mu1=0.01, mu2=0.15):
        ch = self.ch
        return self.da.three_intensity_decoy(
            mu1, mu2, ch.vacuum_yield(),
            ch.intensity_yield(mu1, 500), ch.intensity_yield(mu2, 500),
            ch.qber_x_basis(mu1, 500))
    def test_s1_in_unit_interval(self):    assert 0.0 <= self._decoy().s1 <= 1.0
    def test_e1ph_in_valid_range(self):    assert 0.0 <= self._decoy().e1_ph <= 0.5
    def test_validate_bounds_passes(self): assert self.da.validate_bounds(self._decoy())
    def test_higher_eta_gives_higher_s1(self):
        def s1(eta):
            ch = FakeChannel(eta=eta)
            return self.da.three_intensity_decoy(
                0.01, 0.15, ch.vacuum_yield(),
                ch.intensity_yield(0.01, 500), ch.intensity_yield(0.15, 500),
                ch.qber_x_basis(0.01, 500)).s1
        assert s1(1e-2) >= s1(1e-4)
    def test_mu1_lt_mu2_constraint(self):
        r = self._decoy(mu1=0.15, mu2=0.01)
        assert r.s1 >= 0.0


class TestBinaryEntropy:
    def test_boundary_zero(self):    assert binary_entropy(0.0) == pytest.approx(0.0)
    def test_boundary_one(self):     assert binary_entropy(1.0) == pytest.approx(0.0)
    def test_maximum_at_half(self):  assert binary_entropy(0.5) == pytest.approx(1.0, rel=1e-6)
    def test_symmetry(self):
        for x in [0.1, 0.2, 0.3, 0.4]:
            assert binary_entropy(x) == pytest.approx(binary_entropy(1.0 - x), rel=1e-10)
    def test_monotone_increasing_to_half(self):
        xs = np.linspace(0.01, 0.5, 50)
        hs = [binary_entropy(x) for x in xs]
        assert all(hs[i] <= hs[i+1] for i in range(len(hs)-1))
    def test_clip_out_of_range(self):
        assert binary_entropy(-0.1) == pytest.approx(0.0)
        assert binary_entropy(1.1)  == pytest.approx(0.0)


class TestKeyRate:
    def setup_method(self): self.calc = SNSKeyRateCalculator(f=1.1)
    def test_zero_s1_gives_zero_rate(self):
        assert self.calc.key_rate(0.05, 0.05, 0.0, 0.1, 1e-5, 0.0) == pytest.approx(0.0)
    def test_high_phase_error_gives_zero_rate(self):
        assert self.calc.key_rate(0.05, 0.05, 1e-4, 0.5, 1e-5, 0.0) == 0.0
    def test_rate_non_negative(self):
        for eps in [0.01, 0.1, 0.3]:
            assert self.calc.key_rate(eps, 0.05, 1e-4, 0.05, 1e-6, 0.0) >= 0.0
    def test_rate_decreases_with_e1ph(self):
        r1 = self.calc.key_rate(0.1, 0.1, 1e-3, 0.05, 1e-5, 0.0)
        r2 = self.calc.key_rate(0.1, 0.1, 1e-3, 0.20, 1e-5, 0.0)
        assert r1 >= r2
    def test_key_rate_vs_distance_returns_list(self):
        results = self.calc.key_rate_vs_distance(
            np.linspace(0, 500, 10), FakeChannel(eta=1e-3), 0.05, 0.05, 0.01, 0.15)
        assert len(results) == 10
    def test_key_rate_decreases_with_distance(self):
        ch = FakeChannel(eta=1e-3)
        rates = [r.key_rate for r in self.calc.key_rate_vs_distance(
            np.linspace(10, 300, 20), ch, 0.05, 0.05, 0.01, 0.15)]
        assert all(rates[i] >= rates[i+1] - 1e-15 for i in range(len(rates)-1))


class TestIntegration:
    def test_full_pipeline_runs(self):
        ch = FakeChannel(eta=1e-3)
        da = DecoyStateAnalyzer()
        decoy = da.three_intensity_decoy(
            0.01, 0.15, ch.vacuum_yield(),
            ch.intensity_yield(0.01, 500), ch.intensity_yield(0.15, 500),
            ch.qber_x_basis(0.01, 500))
        assert da.validate_bounds(decoy)
        rate = SNSKeyRateCalculator().key_rate(
            0.05, 0.05, decoy.s1, decoy.e1_ph, ch.intensity_yield(0.05, 500))
        assert rate >= 0.0

    def test_higher_misalignment_reduces_rate(self):
        da   = DecoyStateAnalyzer()
        calc = SNSKeyRateCalculator()
        def get_rate(ea):
            ch = FakeChannel(eta=1e-3, ea=ea)
            decoy = da.three_intensity_decoy(
                0.01, 0.15, ch.vacuum_yield(),
                ch.intensity_yield(0.01, 500), ch.intensity_yield(0.15, 500),
                ch.qber_x_basis(0.01, 500))
            return calc.key_rate(0.05, 0.05, decoy.s1, decoy.e1_ph,
                                 ch.intensity_yield(0.05, 500))
        assert get_rate(0.0) >= get_rate(0.25)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])