"""Slant-path turbulence & scintillation module (notebook Part V / Phase 4).

Hufnagel-Valley Cn2(h) vertical profile, Fried parameter and (corrected,
spherical-wave/uplink) Rytov variance as path integrals through that profile,
plus provisional angle-of-arrival jitter and Strehl-ratio estimates for the
Phase 5 PAT budget. Reuses (does not rebuild) the Gamma-Gamma machinery of
Part II.

Origin: notebook cells 116 (cn2_hv, LAMBDA_PHASE4_M/K_810/H_TOP_M), 118
(fried_r0), 120 (SPHERICAL_UPLINK_FACTOR, sigma_R2_slant_planewave,
sigma_R2_slant -- "Priority 1" uplink correction), 122
(TURBULENCE_SITE_CASES), 128 (D_PROVISIONAL_M, aoa_jitter_urad,
strehl_uncompensated).
"""

import numpy as np
from scipy.integrate import quad

LAMBDA_PHASE4_M = 810e-9   # match Phase 2/3 operating wavelength
K_810 = 2*np.pi/LAMBDA_PHASE4_M
H_TOP_M = 30000.0          # turbulence integral truncation altitude (HV profile negligible above this)

# Andrews Field Guide horizontal spherical/plane ratio: spherical-wave horizontal
# 0.5*Cn2*k^(7/6)*L^(11/6); plane-wave horizontal 1.23*Cn2*k^(7/6)*L^(11/6).
# Matches Doolittle et al. arXiv:2412.03356 Eq. B7 (beta_0^2 = 0.4065 sigma_R^2).
SPHERICAL_UPLINK_FACTOR = 0.5 / 1.23   # = 0.4065

D_PROVISIONAL_M = 1.0   # Micius ground receiver, placeholder pending Phase 5

TURBULENCE_SITE_CASES = {
    "shared baseline (HV 5/7)":            dict(A0=1.7e-14, v=21.0, confidence="Medium (literature HV5/7 default; no city-specific measurement located)"),
    "Delhi urban heat-island (indicative)": dict(A0=3.4e-14, v=21.0, confidence="Indicative (illustrative x2 sensitivity case, not a measured value)"),
}


def cn2_hv(h_m, A0=1.7e-14, v=21.0):
    """Hufnagel-Valley Cn2(h) profile, h in metres, A0 in m^(-2/3), v in m/s
    (rms wind). Cell 116."""
    h = np.asarray(h_m, dtype=float)
    term1 = 0.00594*(v/27.0)**2 * (h*1e-5)**10 * np.exp(-h/1000.0)
    term2 = 2.7e-16*np.exp(-h/1500.0)
    term3 = A0*np.exp(-h/100.0)
    return term1 + term2 + term3


def fried_r0(zenith_deg, lam_m, A0=1.7e-14, v=21.0, h_top=H_TOP_M):
    """Fried parameter via the Fried (1966) path-integral formula. Cell 118."""
    zenith_deg = np.clip(zenith_deg, 0, 89.9)   # below-horizon values are unphysical for this model
    k = 2*np.pi/lam_m
    sec_zeta = 1.0/np.cos(np.radians(zenith_deg))
    integral, _ = quad(lambda h: cn2_hv(h, A0, v), 0, h_top, limit=400)
    return (0.423 * k**2 * sec_zeta * integral)**(-3/5), integral


def sigma_R2_slant_planewave(zenith_deg, lam_m, cn2_profile_fn, h0=0.0, h_top=H_TOP_M):
    """General downlink (plane-wave) path-integral Rytov variance through an
    arbitrary Cn2(h) profile. Cell 120."""
    zenith_deg = np.clip(zenith_deg, 0, 89.9)   # below-horizon values are unphysical for this model
    k = 2*np.pi/lam_m
    sec_zeta = 1.0/np.cos(np.radians(zenith_deg))
    integral, _ = quad(lambda h: cn2_profile_fn(h)*(h-h0)**(5/6), h0, h_top, limit=400)
    return 2.25 * k**(7/6) * sec_zeta**(11/6) * integral


def sigma_R2_slant(zenith_deg, lam_m, cn2_profile_fn, h0=0.0, h_top=H_TOP_M):
    """CORRECTED (Priority 1): spherical-wave (uplink) Rytov variance. SNS-TF-QKD
    transmits ground -> satellite (uplink), so this is the quantity feeding
    gamma_gamma_params() throughout Parts V-XIII, not the plane-wave value.
    Cell 120."""
    return SPHERICAL_UPLINK_FACTOR * sigma_R2_slant_planewave(zenith_deg, lam_m, cn2_profile_fn, h0, h_top)


def aoa_jitter_urad(zenith_deg, D=D_PROVISIONAL_M, lam_m=LAMBDA_PHASE4_M):
    """Angle-of-arrival jitter (Tyler 1994 tracking-error formula). Cell 128."""
    zenith_deg = np.clip(zenith_deg, 0, 89.9)
    sec_zeta = 1.0/np.cos(np.radians(zenith_deg))
    integral, _ = quad(lambda h: cn2_hv(h), 0, H_TOP_M, limit=400)
    sigma_theta2 = 2.914 * D**(-1/3) * sec_zeta * integral
    return np.sqrt(sigma_theta2)*1e6


def strehl_uncompensated(r0, D=D_PROVISIONAL_M):
    """Marechal/Noll approximation, D/r0 scaling. Cell 128."""
    return np.exp(-1.03*(D/r0)**(5/3))
