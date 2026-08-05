"""Generalized pass construction and per-link real-loss/joint-optimized
key-rate evaluation (notebook Part XII / Phase 9).

Reuses every already-validated function from Parts III-X unchanged; only the
per-city input parameters differ. Scoped to the geometry + atmosphere +
turbulence + real-loss + joint-optimized key-rate depth of Parts III-VI/X --
the full seasonal/weather-gated annual analysis is not repeated here (see
notebook Section 12.7).

Origin: notebook cell 229 (build_pass), cell 231 (evaluate_link).
"""

import numpy as np

from src.satellite.geometry import (
    OMEGA_E, geodetic_to_ecef, gc_midpoint, eci_position, eci_to_ecef, look_angles, haversine_km,
)
from src.atmosphere.transmittance import LAMBDA_UM, TAU_GAS_ZENITH, kasten_young_airmass, tau_aerosol_zenith, tau_rayleigh_zenith
from src.atmosphere.turbulence import sigma_R2_slant, cn2_hv, LAMBDA_PHASE4_M
from src.atmosphere.seasonal import A0_BASELINE
from src.pat.losses import (
    ETA_RX_INTERNAL, ETA_DET_REAL, E_A_MISALIGNMENT, R_TX_HZ, R_NOISE_BASELINE_HZ, P_DARK_REAL,
    eta_diffraction,
)
from src.optimization.joint import optimize_joint_ensemble
from src.links.city_pairs import CITY_DB


def build_pass(city_a, city_b, orb, a, inc, gmst0, el_min=20.0, t_half=500.0, n_t=2001):
    """Cell 229: engineer a great-circle-midpoint crossing for any two
    cities, using the fixed orbital plane from Part III's ``construct_pass``
    (``orb``, ``a``, ``inc``, ``gmst0``)."""
    a_p, b_p = CITY_DB[city_a], CITY_DB[city_b]
    dist_km = haversine_km(a_p['lat'], a_p['lon'], b_p['lat'], b_p['lon'])

    ecef_a = geodetic_to_ecef(a_p['lat'], a_p['lon'], a_p['alt_km'])
    ecef_b = geodetic_to_ecef(b_p['lat'], b_p['lon'], b_p['alt_km'])

    lat_m, lon_m_ = gc_midpoint(a_p['lat'], a_p['lon'], b_p['lat'], b_p['lon'])
    lat_t, lon_t = np.radians(lat_m), np.radians(lon_m_)
    u0_ = np.arcsin(np.clip(np.sin(lat_t)/np.sin(inc), -1, 1))
    lon_arg = np.arctan2(np.cos(inc)*np.sin(u0_), np.cos(u0_))
    raan0_ = lon_t + gmst0 - lon_arg

    t_grid_ = np.linspace(-t_half, t_half, n_t)
    el_a = np.empty_like(t_grid_); el_b = np.empty_like(t_grid_)
    rng_a = np.empty_like(t_grid_); rng_b = np.empty_like(t_grid_)
    for k, t in enumerate(t_grid_):
        raan_t = raan0_ + orb['raan_dot']*t
        u_t = u0_ + orb['u_dot']*t
        r_eci = eci_position(a, inc, raan_t, u_t)
        gmst_t = gmst0 + OMEGA_E*t
        r_ecef = eci_to_ecef(r_eci, gmst_t)
        e1, _, rg1 = look_angles(r_ecef, ecef_a, a_p['lat'], a_p['lon'])
        e2, _, rg2 = look_angles(r_ecef, ecef_b, b_p['lat'], b_p['lon'])
        el_a[k], el_b[k] = e1, e2
        rng_a[k], rng_b[k] = rg1, rg2

    common = (el_a >= el_min) & (el_b >= el_min)
    return dict(dist_km=dist_km, t=t_grid_, el_a=el_a, el_b=el_b, rng_a=rng_a, rng_b=rng_b,
                common=common, dt=t_grid_[1]-t_grid_[0])


def evaluate_link(city_a, city_b, pr, eta_point_val, suppression_factor=1000.0, e_a=E_A_MISALIGNMENT):
    """Cell 231: fully-dressed real-loss channel + joint-optimized key rate
    at the link's minimum-joint-slant-range point."""
    idx_common = np.where(pr['common'])[0]
    if len(idx_common) == 0:
        return None
    i_best = idx_common[np.argmin(pr['rng_a'][idx_common] + pr['rng_b'][idx_common])]

    results_pair = {}
    for label, city, rng, el in [("A", city_a, pr['rng_a'], pr['el_a']), ("B", city_b, pr['rng_b'], pr['el_b'])]:
        cp = CITY_DB[city]
        zenith = 90.0 - el[i_best]
        airmass = kasten_young_airmass(zenith)
        tau_a = tau_aerosol_zenith(cp['aod500'], cp['alpha'], 500.0)
        tau_ray = tau_rayleigh_zenith(LAMBDA_UM, P_site_hpa=1013.25*np.exp(-cp['alt_km']/8.5))
        tau_zenith_total = tau_ray + tau_a + TAU_GAS_ZENITH
        eta_atm = np.exp(-tau_zenith_total*airmass)

        L_m = rng[i_best]*1000.0
        eta_diff = eta_diffraction(L_m)
        eta_total = eta_diff * eta_atm * eta_point_val * ETA_RX_INTERNAL * ETA_DET_REAL

        A0_city = A0_BASELINE * cp['turb_mult']
        sigma_R2 = sigma_R2_slant(zenith, LAMBDA_PHASE4_M, lambda h, A0=A0_city: cn2_hv(h, A0=A0))

        p_bg_baseline = R_NOISE_BASELINE_HZ / R_TX_HZ
        p_dark_total = P_DARK_REAL + (p_bg_baseline*2.0)/suppression_factor  # Delhi-level (x2) background, uniform conservative assumption

        R_opt, eps_opt, mu2_opt, mu1_opt = optimize_joint_ensemble(eta_total, e_a, sigma_R2, p_dark_total,
                                                                     n_mc=15000, n_grid=10)
        results_pair[label] = dict(city=city, elevation=el[i_best], slant_range_km=rng[i_best],
                                    tau_zenith_total=tau_zenith_total, eta_atm=eta_atm, eta_total=eta_total,
                                    sigma_R2=sigma_R2, R=R_opt, eps=eps_opt, mu2=mu2_opt, mu1=mu1_opt)
    return results_pair
