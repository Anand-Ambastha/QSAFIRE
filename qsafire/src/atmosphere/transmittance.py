"""Static atmospheric transmittance module (notebook Part IV / Phase 3).

Combines Rayleigh (molecular) scattering, aerosol (Mie) extinction, and a
small fixed gaseous-absorption residual by Beer's law, then scales zenith
optical depth to the slant path via the Kasten & Young (1989) airmass
formula, following the standard gamma = alpha_m + alpha_a + beta_m + beta_a
decomposition used in the free-space-QKD literature.

Origin: notebook cells 100 (wavelength), 102 (tau_rayleigh_zenith), 104
(AEROSOL_SITES, tau_aerosol_zenith), 106 (kim_kruse_q,
tau_aerosol_visibility), 108 (TAU_GAS_ZENITH), 110 (kasten_young_airmass,
SITE_TAU_ZENITH, SITE_TAU_ZENITH_WINTER).
"""

import numpy as np

LAMBDA_UM = 0.810     # operating wavelength, micrometers (810 nm -- matches Micius & Behera/Sinha 2024)
LAMBDA_NM = 810.0

TAU_GAS_ZENITH = 0.02   # indicative small residual, 810 nm atmospheric window

AEROSOL_SITES = {
    "Delhi (Alice) -- annual baseline":  dict(aod0=0.78, alpha=0.78, lam0_nm=500.0, confidence="High (NPL New Delhi, 2007-2008 annual mean)"),
    "Delhi (Alice) -- winter haze case": dict(aod0=0.91, alpha=0.70, lam0_nm=500.0, confidence="High (Ganguly et al. 2006, Dec 2004 campaign; alpha not reported, representative value assumed)"),
    "Mumbai (Bob) -- indicative":        dict(aod0=0.35, alpha=0.90, lam0_nm=500.0, confidence="Indicative (regional 0.2-0.5 range, ISRO-GBP campaign + MODIS multi-city qualitative comparison)"),
}


def tau_rayleigh_zenith(lam_um, P_site_hpa=1013.25, P_std_hpa=1013.25):
    """Sea-level-standard Rayleigh (molecular) zenith optical depth, scaled
    linearly to site surface pressure. Cell 102."""
    A, B, C, D = 0.00864, 3.916, 0.074, 0.05
    tau_std = A * lam_um**(-(B + C*lam_um + D/lam_um))
    return tau_std * (P_site_hpa / P_std_hpa)


def tau_aerosol_zenith(aod0, alpha, lam0_nm, lam_nm=LAMBDA_NM):
    """Angstrom-power-law scaling of column AOD to the operating wavelength. Cell 104."""
    return aod0 * (lam_nm/lam0_nm)**(-alpha)


def kim_kruse_q(V_km):
    """Kim/Kruse size-distribution parameter, selected by visibility regime. Cell 106."""
    if V_km > 50: return 1.6
    elif V_km > 6: return 1.3
    elif V_km > 1: return 0.16*V_km + 0.34
    elif V_km > 0.5: return V_km - 0.5
    else: return 0.0


def tau_aerosol_visibility(V_km, lam_nm=LAMBDA_NM, H_aerosol_km=2.0):
    """Visibility-based (Kim/Kruse) aerosol zenith optical depth cross-check. Cell 106."""
    q = kim_kruse_q(V_km)
    sigma_per_km = (3.91/V_km) * (lam_nm/550.0)**(-q)   # km^-1
    return sigma_per_km * H_aerosol_km, q


def kasten_young_airmass(zenith_deg):
    """Kasten & Young (1989) airmass formula. Cell 110."""
    zenith_deg = np.clip(zenith_deg, 0, 89.9)
    return 1.0/(np.cos(np.radians(zenith_deg)) + 0.50572*(96.07995 - zenith_deg)**(-1.6364))


def _aerosol_kwargs(case):
    p = AEROSOL_SITES[case]
    return {k: v for k, v in p.items() if k in ('aod0', 'alpha', 'lam0_nm')}


def site_tau_zenith(tau_R_810):
    """Total baseline zenith optical depth per station (Rayleigh + aerosol +
    gas). Cell 110."""
    return {
        "Delhi (Alice)": tau_R_810 + tau_aerosol_zenith(**_aerosol_kwargs("Delhi (Alice) -- annual baseline")) + TAU_GAS_ZENITH,
        "Mumbai (Bob)":  tau_R_810 + tau_aerosol_zenith(**_aerosol_kwargs("Mumbai (Bob) -- indicative")) + TAU_GAS_ZENITH,
    }


def site_tau_zenith_winter(tau_R_810):
    """Delhi winter-haze sensitivity case zenith optical depth. Cell 110."""
    return {
        "Delhi (Alice) winter-haze": tau_R_810 + tau_aerosol_zenith(**_aerosol_kwargs("Delhi (Alice) -- winter haze case")) + TAU_GAS_ZENITH,
    }
