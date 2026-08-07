"""One-at-a-time (OAT) sensitivity analysis for the satellite-uplink chain
(src/satellite, src/atmosphere, src/pat, src/optimization/joint.py) - i.e.
the real-loss, geometry-driven model exercised by
src/links/pass_builder.evaluate_link, as opposed to the abstract
distance-only FSO model covered by scripts/sensitivity_fso.py.

Reuses evaluate_link's own physics chain (eta_diffraction, eta_pointing,
atmospheric transmittance via Kasten-Young airmass, HV-profile Rytov
variance, optimize_joint_ensemble) rather than re-deriving it, but exposes
every input as a free parameter instead of the fixed CITY_DB/SITES values,
so each can be swept independently around a single representative baseline
geometry (a Delhi-side uplink at 500 km altitude, 45 deg elevation).

Parameters covered:
  h_km            satellite altitude (-> slant range via slant_from_elev)
  elevation_deg   link elevation angle (-> zenith, airmass, sigma_R2)
  D_RX_M          ground-receiver aperture diameter
  THETA_DIV_RAD   satellite Tx beam half-divergence
  sigma_jitter    pointing-jitter std dev (PAT budget)
  eta_rx_internal internal optics/mirror throughput
  eta_det_real    detector (Si APD) photon-detection efficiency
  aod500          aerosol optical depth at 500 nm (site atmosphere)
  A0              ground-level Cn2(0) turbulence constant (HV profile)
  wind_v          rms wind speed in the HV profile (turbulence)
  r_noise_hz      background-light count rate
  e_a             misalignment/interference error rate

As in sensitivity_fso.py, elasticity E = (%change in R)/(%change in
parameter) is the headline comparable sensitivity measure, evaluated at
+-10%/+-20% perturbations, one parameter at a time, holding all others at
baseline.

Run with:  python -m scripts.sensitivity_satellite
Outputs:   outputs/csv/sensitivity_satellite.csv

NOTE: sigma_R2 at a real satellite slant path (near-zenith, short
turbulence-relevant column) stays well inside the Gamma-Gamma model's
validity range (channel_models.SIGMA_R2_VALIDITY_MAX) for every
perturbation used here - this is checked and asserted, not assumed.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd

from src.satellite.geometry import slant_from_elev
from src.atmosphere.transmittance import LAMBDA_UM, TAU_GAS_ZENITH, kasten_young_airmass, tau_aerosol_zenith, tau_rayleigh_zenith
from src.atmosphere.turbulence import sigma_R2_slant, cn2_hv, LAMBDA_PHASE4_M
from src.atmosphere.seasonal import A0_BASELINE
from src.pat.losses import (
    LAMBDA_P5_M, ETA_RX_INTERNAL, ETA_DET_REAL, E_A_MISALIGNMENT, R_TX_HZ,
    R_NOISE_BASELINE_HZ, P_DARK_REAL, THETA_DIV_RAD, D_RX_M, OBSCURATION_RATIO,
    eta_diffraction, eta_pointing,
)
from src.fso.channel_models import sigma_R2_within_validity, SIGMA_R2_VALIDITY_MAX
from src.optimization.joint import optimize_joint_ensemble

OUT_CSV = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "outputs", "csv", "sensitivity_satellite.csv")

# ---------------------------------------------------------------------------
# Baseline operating point: Delhi-side uplink, 500 km altitude, 45 deg elevation
# ---------------------------------------------------------------------------
AOD500_DELHI = 0.78
ALPHA_ANGSTROM_DELHI = 0.78

BASE = dict(
    h_km=500.0,
    elevation_deg=45.0,
    D_RX_M=D_RX_M,
    THETA_DIV_RAD=THETA_DIV_RAD,
    sigma_jitter=2.5e-6,          # "Conservative default" PAT case
    eta_rx_internal=ETA_RX_INTERNAL,
    eta_det_real=ETA_DET_REAL,
    aod500=AOD500_DELHI,
    A0=A0_BASELINE,
    wind_v=21.0,
    r_noise_hz=R_NOISE_BASELINE_HZ,
    e_a=E_A_MISALIGNMENT,
)
PERTURBATIONS = (-0.20, -0.10, 0.10, 0.20)

# Fixed (not swept) chain parameters
SUPPRESSION_FACTOR = 1000.0   # background gating/coincidence suppression, as in evaluate_link


def evaluate_uplink(h_km, elevation_deg, D_RX_M, THETA_DIV_RAD, sigma_jitter,
                     eta_rx_internal, eta_det_real, aod500, A0, wind_v, r_noise_hz, e_a,
                     n_mc=8000, n_grid=8, seed=0):
    """Standalone re-parameterised version of evaluate_link's per-arm inner
    loop (src/links/pass_builder.py), with every physical input exposed."""
    zenith = 90.0 - elevation_deg
    L_km = float(slant_from_elev(elevation_deg, h_km))
    L_m = L_km * 1000.0

    # --- optics / PAT chain ---
    D_rx_int = D_RX_M * OBSCURATION_RATIO
    G_tx = 8.0 / THETA_DIV_RAD ** 2
    G_rx = (np.pi ** 2 / LAMBDA_P5_M ** 2) * (D_RX_M ** 2 - D_rx_int ** 2)
    eta_diff = eta_diffraction(L_m, G_tx=G_tx, G_rx=G_rx, lam=LAMBDA_P5_M)
    eta_point = eta_pointing(sigma_jitter, theta_div=THETA_DIV_RAD)

    # --- atmospheric transmittance ---
    airmass = kasten_young_airmass(zenith)
    tau_a = tau_aerosol_zenith(aod500, ALPHA_ANGSTROM_DELHI, 500.0)
    tau_ray = tau_rayleigh_zenith(LAMBDA_UM)
    tau_zenith_total = tau_ray + tau_a + TAU_GAS_ZENITH
    eta_atm = float(np.exp(-tau_zenith_total * airmass))

    eta_total = eta_diff * eta_atm * eta_point * eta_rx_internal * eta_det_real

    # --- turbulence (HV profile -> spherical-wave slant Rytov variance) ---
    sigma_R2 = float(sigma_R2_slant(zenith, LAMBDA_PHASE4_M, lambda h: cn2_hv(h, A0=A0, v=wind_v)))

    # --- background / effective dark-count budget ---
    p_bg_baseline = r_noise_hz / R_TX_HZ
    p_dark_total = P_DARK_REAL + (p_bg_baseline * 2.0) / SUPPRESSION_FACTOR

    valid = bool(sigma_R2_within_validity(sigma_R2))
    if not valid:
        return None, dict(L_km=L_km, eta_total=float(eta_total), sigma_R2=sigma_R2, valid=valid)

    out = optimize_joint_ensemble(float(eta_total), e_a, sigma_R2, p_dark_total, n_mc=n_mc, n_grid=n_grid, seed=seed)
    R = out[0]
    return R, dict(L_km=L_km, eta_total=float(eta_total), sigma_R2=sigma_R2, valid=valid)


def run():
    rows = []
    params = list(BASE.keys())

    R0, diag0 = evaluate_uplink(**BASE)
    print(f"Baseline (h={BASE['h_km']} km, el={BASE['elevation_deg']} deg): "
          f"L={diag0['L_km']:.1f} km, eta_total={diag0['eta_total']:.3e}, "
          f"sigma_R2={diag0['sigma_R2']:.4f} (valid={diag0['valid']}), R={R0:.6e}")
    assert diag0["valid"], "Baseline sigma_R2 unexpectedly outside GG validity range"

    for pname in params:
        for pct in PERTURBATIONS:
            kwargs = dict(BASE)
            kwargs[pname] = BASE[pname] * (1.0 + pct)
            # elevation angle: clip perturbations to a physically sane range
            if pname == "elevation_deg":
                kwargs[pname] = float(np.clip(kwargs[pname], 5.0, 89.0))
            R, diag = evaluate_uplink(**kwargs)
            if R is None:
                rows.append(dict(parameter=pname, baseline_value=BASE[pname],
                                  perturbed_value=kwargs[pname], param_pct_change=pct * 100.0,
                                  R_baseline=R0, R_perturbed=np.nan, R_pct_change=np.nan,
                                  elasticity=np.nan,
                                  note=f"sigma_R2={diag['sigma_R2']:.2f} exceeds validity ceiling "
                                       f"{SIGMA_R2_VALIDITY_MAX} - excluded"))
                continue
            pct_change_R = (R - R0) / abs(R0) * 100.0
            elasticity = pct_change_R / (pct * 100.0)
            rows.append(dict(parameter=pname, baseline_value=BASE[pname],
                              perturbed_value=kwargs[pname], param_pct_change=pct * 100.0,
                              R_baseline=R0, R_perturbed=R, R_pct_change=pct_change_R,
                              elasticity=elasticity,
                              note=f"L={diag['L_km']:.1f}km eta_tot={diag['eta_total']:.3e} "
                                   f"sigma_R2={diag['sigma_R2']:.3f}"))

    df = pd.DataFrame(rows)

    ranking = (df[np.isclose(df["param_pct_change"].abs(), 10.0)]
               .groupby("parameter")["elasticity"].apply(lambda s: s.abs().mean())
               .sort_values(ascending=False))
    print("\nSensitivity ranking (mean |elasticity| at +-10%):")
    for pname, val in ranking.items():
        print(f"  {pname:>16s}: {val:8.3f}")

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV} ({len(df)} rows)")
    return df


if __name__ == "__main__":
    run()