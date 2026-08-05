"""Smoke tests for the migrated Part I (fibre) and Part II (FSO/turbulence)
modules.

These are lightweight numerical-regression checks against values printed by
the original notebook (see notebooks/original_notebook.ipynb), not a full
test suite. Run with:

    python -m tests.smoke_test
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from src.simulation.fiber_channel import channel_eta, sns_key_rate, optimize_key_rate
from src.fso.channel_models import fso_deterministic_eta, gamma_gamma_params, rytov_variance
from src.fso.quadrature import gauss_laguerre_ensemble_rate
from src.fso.monte_carlo import monte_carlo_ensemble_rate
from src.fso.config import CN2_MODERATE, E_A
from src.satellite.geometry import slant_from_elev, sun_sync_orbit
from src.atmosphere.transmittance import tau_rayleigh_zenith
from src.atmosphere.turbulence import fried_r0
from src.links.city_pairs import CITY_DB
from src.links.metrics import validate_symmetric_limit


def test_channel_eta_matches_notebook():
    # Section 3 sanity table: eta_TFQKD(L=100) should equal 0.08 (printed in
    # the notebook's Section 3 illustration).
    assert np.isclose(channel_eta(100), 0.08, rtol=1e-6)
    print("PASS: channel_eta(100 km) == 0.08")


def test_fiber_key_rate_quick_check():
    # Notebook Section 8 quick check:
    # R(L=200km, ea=0.15, eps=0.05, mu'=0.3) = -2.2872139670263216e-05
    R = sns_key_rate(0.05, 0.3, 200, 0.15)
    assert np.isclose(R, -2.2872139670263216e-05, rtol=1e-6)
    print(f"PASS: sns_key_rate(0.05, 0.3, 200, 0.15) = {R:.6e}")


def test_optimizer_finds_positive_rate_at_300km_25pct():
    R_opt, eps_opt, mu_opt = optimize_key_rate(300, 0.25)
    assert R_opt > 0
    print(f"PASS: optimize_key_rate(300, 0.25) -> R={R_opt:.4e} > 0")


def test_fso_deterministic_eta_matches_fiber_channel_eta():
    # Notebook Part II intro: fso_deterministic_eta(d_km) with the FIBRE alpha
    # should equal channel_eta(2*d_km) - both express eta_d * 10^(-alpha*d_arm/10)
    # with the SAME per-arm distance convention.
    d = 150.0
    eta_fso_fiber_alpha = fso_deterministic_eta(d, alpha=0.2, eta_d=0.8)
    eta_fiber = channel_eta(2 * d, alpha=0.2, eta_d=0.8)
    assert np.isclose(eta_fso_fiber_alpha, eta_fiber, rtol=1e-12)
    print("PASS: fso_deterministic_eta(d, fibre-alpha) == channel_eta(2d)")


def test_gamma_gamma_params_unit_mean():
    ag, bg = gamma_gamma_params(3.0)
    assert ag > 0 and bg > 0
    print(f"PASS: gamma_gamma_params(3.0) -> alpha_g={ag:.4f}, beta_g={bg:.4f}")


def test_quadrature_matches_monte_carlo():
    d = 6.0
    s2 = rytov_variance(d, CN2_MODERATE)
    R_gl = gauss_laguerre_ensemble_rate(0.02, 0.3, d, E_A, s2, n_points=40)
    mc = monte_carlo_ensemble_rate(0.02, 0.3, d, E_A, s2, n_mc=50_000,
                                    rng=np.random.default_rng(0))
    rel_err = abs(R_gl - mc["mean"]) / abs(mc["mean"])
    assert rel_err < 0.02, f"relative error {rel_err} exceeds smoke-test tolerance"
    print(f"PASS: Gauss-Laguerre vs Monte Carlo relative error = {rel_err:.4%}")


def test_slant_from_elev_matches_liao_2017():
    # Notebook Part III, cell 89: elevation<->slant-range formula validated
    # against Liao et al. (2017), within 2%.
    d1 = slant_from_elev(85.7, 500.0)
    d2 = slant_from_elev(25.0, 500.0)
    assert abs(100*(d1-507.0)/507.0) < 2.0
    assert abs(100*(d2-1034.7)/1034.7) < 2.0
    print(f"PASS: slant_from_elev matches Liao et al. (2017) within 2% ({d1:.1f} km, {d2:.1f} km)")


def test_sun_sync_orbit_solves_landsat_inclination():
    # Landsat-8/9 (705 km) real inclination is 98.2 deg.
    orb = sun_sync_orbit(705.0)
    inc_deg = np.degrees(orb['i'])
    assert abs(inc_deg - 98.2) < 0.5
    print(f"PASS: sun_sync_orbit(705 km) inclination = {inc_deg:.2f} deg (Landsat-8/9 real value: 98.2 deg)")


def test_rayleigh_matches_ams_glossary():
    # Notebook Part IV, cell 102: Rayleigh zenith optical depth @ 500 nm vs.
    # AMS Glossary of Meteorology reference value 0.145.
    tau = tau_rayleigh_zenith(0.500)
    assert abs(100*(tau-0.145)/0.145) < 5
    print(f"PASS: tau_rayleigh_zenith(500nm) = {tau:.4f} (AMS Glossary reference: 0.145)")


def test_fried_r0_matches_hv57_canonical():
    # Notebook Part V, cell 118: HV 5/7 canonical zenith Fried parameter target 5.00 cm.
    r0, _ = fried_r0(0.0, 500e-9)
    assert abs(100*(r0*100-5.0)/5.0) < 5
    print(f"PASS: fried_r0(zenith, 500nm) = {r0*100:.2f} cm (HV 5/7 canonical target: 5.00 cm)")


def test_asymmetric_rate_formulas_reduce_to_symmetric():
    # Notebook Part XIII, cell 240: eta_A=eta_B must reproduce the original
    # symmetric formulas to machine precision.
    max_dev = validate_symmetric_limit()
    assert max_dev < 1e-12
    print(f"PASS: asymmetric rate formulas reduce to symmetric case (max deviation {max_dev:.2e})")


def test_city_db_has_six_cities():
    assert len(CITY_DB) == 6
    print(f"PASS: CITY_DB has {len(CITY_DB)} cities")


def run_all():
    tests = [v for k, v in globals().items() if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
    print(f"\nAll {len(tests)} smoke tests passed.")


if __name__ == "__main__":
    run_all()
