# API Reference

Function-level reference for every public function in the package, grouped
by package. All angles are in **degrees** unless a parameter name ends in
`_rad`; all distances are in **km** unless a parameter name ends in `_m`;
all wavelengths follow the module's own convention (noted per module). For
the physics behind these formulas, see `docs/theory.md`.

---

## `src/engine.py` — canonical SNS-TF-QKD rate equations

The single source of truth for the protocol's rate formulas. Every other
channel model (fiber, FSO, satellite, multi-city links) calls into this
module (or, for asymmetric links, `src/links/metrics.py`'s direct
generalization of it) rather than re-deriving the protocol.

Module constants: `P_DARK` (default per-pulse dark-count probability),
`F_EC` (error-correction inefficiency factor), `DECOY_RATIO` (default
mu1/mu2 ratio for the fixed-ratio formula).

| Function | Signature | Returns |
|---|---|---|
| `p0_weight` | `(mu)` | Vacuum-state Poisson weight at intensity `mu`. |
| `p1_weight` | `(mu)` | Single-photon Poisson weight at intensity `mu`. |
| `p2_weight` | `(mu)` | Two-photon Poisson weight at intensity `mu`. |
| `s0_vacuum_yield` | `(p_dark=P_DARK)` | Vacuum-state yield `s0` from the dark-count probability. |
| `gain_X` | `(mu, eta, p_dark=P_DARK)` | X-window gain `S_mu` at intensity `mu`, transmittance `eta`. |
| `error_X` | `(mu, eta, e_a, p_dark=P_DARK)` | X-window QBER given misalignment error `e_a`. |
| `S_Z_gain` | `(eps, mu_prime, eta, p_dark=P_DARK)` | Z-window gain, sending probability `eps`. |
| `E_Z_qber` | `(eps, mu_prime, eta, p_dark=P_DARK)` | Z-window QBER (dark-count dominated). |
| `s1_decoy_lower_bound` | `(mu1, mu2, S_mu1, S_mu2, s0)` | Decoy-state lower bound on the single-photon yield `s1`. |
| `e1ph_upper_bound` | `(mu1, S_mu1, E_mu1, s0, s1)` | Decoy-state upper bound on the phase-error rate `e1ph`. |
| `binary_entropy` | `(x)` | Shannon binary entropy `h2(x)`, safe at `x=0,1`. |
| `sns_key_rate_eta` | `(eps, mu_prime, eta, e_a, decoy_ratio=DECOY_RATIO, f_ec=F_EC, p_dark=P_DARK)` | Asymptotic secure key rate (bits/pulse), `eta`-driven, fixed decoy ratio (`mu1 = decoy_ratio*mu2`). |
| `full_metrics_eta_vec` | same args | Vectorized version returning a dict of every intermediate quantity (`S1,S2,E1,s0,s1,e1ph,Sz,Ez,R`) for array `eta`. |
| `full_metrics_eta` | same args | Scalar convenience wrapper around `full_metrics_eta_vec`. |

---

## `src/simulation/fiber_channel.py` — fiber channel

| Function | Signature | Returns |
|---|---|---|
| `channel_eta` | `(L_km, alpha=ALPHA, eta_d=ETA_D)` | Fiber transmittance at distance `L_km`: `eta_d * 10**(-alpha*L_km/10)`. |
| `true_single_photon_yield` | `(eta, p_dark=P_DARK)` | Reference (non-decoy-bounded) single-photon yield, for validation. |
| `sns_key_rate` | `(eps, mu_prime, L, e_a, decoy_ratio=0.2, f_ec=F_EC, p_dark=P_DARK)` | Key rate as a function of fiber distance `L` (wraps `channel_eta` + `src.engine.sns_key_rate_eta`). |
| `optimize_key_rate` | `(L, e_a, decoy_ratio=0.2, ...)` | Grid + local optimization over `(eps, mu_prime)` at fixed `L`; returns `(R_opt, eps_opt, mu_opt)`. |

---

## `src/fso/` — free-space-optical channel with Gamma-Gamma turbulence

### `channel_models.py`

| Function | Signature | Returns |
|---|---|---|
| `fso_deterministic_eta` | `(d_km, alpha=ALPHA_FSO_DB_KM, eta_d=0.8)` | Deterministic (no-turbulence) FSO transmittance. |
| `channel_loss_dB` | `(d_km, alpha=ALPHA_FSO_DB_KM)` | Deterministic loss in dB. |
| `rytov_variance` | `(d_km, Cn2, k=K_WAVENUMBER)` | Horizontal-path (plane-wave) Rytov variance `sigma_R^2` for constant `Cn2`. |
| `gamma_gamma_params` | `(sigma_R2)` | `(alpha_g, beta_g)` Gamma-Gamma shape parameters from `sigma_R2`. |
| `classify_turbulence_regime` | `(sigma_R2)` | `"weak"` / `"moderate"` / `"strong"` string (array-safe) from `sigma_R2`. |
| `validity_flag` | `(sigma_R2)` | Whether the weak-fluctuation Gamma-Gamma model is expected to hold. |
| `gamma_gamma_pdf` | `(h, alpha_g, beta_g)` | Gamma-Gamma irradiance PDF value at `h`. |
| `sample_gamma_gamma` | `(alpha_g, beta_g, size, rng)` | Draws `size` samples of unit-mean Gamma-Gamma fading (product of two Gamma variates). |

### `monte_carlo.py`

| Function | Signature | Returns |
|---|---|---|
| `monte_carlo_ensemble_rate` | `(eps, mu_prime, d_km, e_a, sigma_R2, n_mc=N_MC_DEFAULT, ...)` | Mean key rate under turbulence fading, by direct Monte Carlo sampling. |
| `monte_carlo_convergence` | `(eps, mu_prime, d_km, e_a, sigma_R2, n_list, ...)` | Mean/std of the Monte Carlo estimate across a list of sample sizes, for a convergence study. |

### `quadrature.py`

| Function | Signature | Returns |
|---|---|---|
| `gauss_laguerre_ensemble_rate` | `(eps, mu_prime, d_km, e_a, sigma_R2, n_points=40, ...)` | Mean key rate under turbulence fading, by Gauss-Laguerre quadrature (closed-form alternative to Monte Carlo). |
| `gauss_laguerre_ensemble_metrics` | same args | Full metrics dict version of the above. |
| `quadrature_validation` | `(eps, mu_prime, d_km, e_a, sigma_R2, n_mc=N_MC_DEFAULT, ...)` | Cross-checks the quadrature result against Monte Carlo; returns both plus the relative error. |

### `optimization.py`

| Function | Signature | Returns |
|---|---|---|
| `optimize_instantaneous` | `(eta_value, e_a, eps_bounds=EPS_BOUNDS, mu_bounds=MU_BOUNDS, ...)` | Grid + local `(eps, mu')` optimization at one deterministic `eta` value. |
| `instantaneous_maps` | `(distances_km, e_a, sigma_R2_fn, ...)` | Builds instantaneous-optimum maps (rate, eps*, mu*) across a distance grid. |
| `grid_search_joint_ensemble` | `(d_km, e_a, sigma_R2, n_grid=20, n_points=20, ...)` | Coarse grid search over `(eps, mu')` under ensemble (turbulence-faded) averaging. |
| `differential_evolution_joint_ensemble` | `(d_km, e_a, sigma_R2, n_points=20, ...)` | Global optimizer (SciPy Differential Evolution) over the ensemble rate. |
| `basin_hopping_joint_ensemble` | `(d_km, e_a, sigma_R2, x0, n_points=20, ...)` | Global optimizer (SciPy Basin Hopping) seeded at `x0`. |
| `joint_optimize_ensemble` | `(d_km, e_a, sigma_R2, n_points=20, ...)` | Triangulates grid search + DE + basin hopping and returns the best result found. |
| `optimize_epsilon_ensemble` | `(d_km, mu_fixed, e_a, sigma_R2, n_points=40, ...)` | 1-D optimization over `eps` only, at fixed `mu'`. |

### `metrics.py`

| Function | Signature | Returns |
|---|---|---|
| `build_deterministic_table` | `(distances_km, e_a=E_A)` | DataFrame of the deterministic-channel key rate/QBER/etc. across a distance grid. |
| `build_gamma_gamma_table` | `(distances_km, turbulence_cases=TURBULENCE_CASES, e_a=E_A, n_mc=100_000, seed=0)` | DataFrame of joint-optimized ensemble results across distances x turbulence cases. |
| `maximum_performance_summary` | `(df_det, df_gg)` | Summary dict/DataFrame of max SKR, max secure distance, max tolerable loss per channel model. |

### `config.py`, `visualization.py`, `export.py`

`config.py` holds shared constants only (wavelength, attenuation, distance
grids, turbulence-case definitions, optimization bounds, output paths — no
functions). `visualization.py` and `export.py` hold plotting and
CSV/report-writing routines called from `src/experiments/fso_experiments.py`
rather than used standalone.

---

## `src/satellite/geometry.py` — geodesy, orbit design, propagation

Module constants: `MU_EARTH`, `RE_EQ`, `F_WGS84`, `E2_WGS84`, `J2`,
`OMEGA_E`, `YEAR_TROPICAL_DAYS` (standard geodetic/gravitational constants);
`SITES` (the Delhi/Mumbai ground-station dict).

| Function | Signature | Returns |
|---|---|---|
| `geodetic_to_ecef` | `(lat_deg, lon_deg, alt_km)` | ECEF position vector (km), WGS-84. |
| `enu_rotation` | `(lat_deg, lon_deg)` | 3x3 matrix whose rows are the East/North/Up unit vectors in ECEF. |
| `build_site_ecef` | `()` | Dict of ECEF vectors for every entry in `SITES`. |
| `sun_sync_orbit` | `(h_km, Re=RE_EQ, mu=MU_EARTH, j2=J2)` | Solves for the sun-synchronous inclination at altitude `h_km`; returns a dict with `a, i, n0, n_true, raan_dot, u_dot`. |
| `julian_date` | `(year, month, day, hour=0, minute=0, sec=0.0)` | Julian date. |
| `gmst_rad` | `(jd)` | Greenwich Mean Sidereal Time (rad) from a Julian date. |
| `eci_position` | `(a, i, raan, u)` | Satellite position in ECI, given semi-major axis, inclination, RAAN, argument of latitude. |
| `eci_to_ecef` | `(r_eci, gmst)` | Rotates an ECI vector into ECEF at a given GMST. |
| `look_angles` | `(sat_ecef, site_ecef_vec, site_lat, site_lon)` | `(elevation_deg, azimuth_deg, slant_range_km)` from a ground site to a satellite position. |
| `refraction_bump_deg` | `(true_elev_deg)` | Bennett (1982) approximate atmospheric refraction correction, degrees. |
| `slant_from_elev` | `(elev_deg, h_km, Re=RE_EQ)` | Spherical-Earth slant range from elevation angle and satellite altitude. |
| `gc_midpoint` | `(lat1, lon1, lat2, lon2)` | Great-circle midpoint of two lat/lon points. |
| `sub_satellite` | `(t, a, inc, raan0, u0, raan_dot, u_dot, gmst0)` | Sub-satellite ground track `(lat_deg, lon_deg)` at time `t` (array-safe). |
| `haversine_km` | `(lat1, lon1, lat2, lon2, R=RE_EQ)` | Great-circle distance. |

---

## `src/atmosphere/transmittance.py` — static zenith optical depth

Module constants: `LAMBDA_UM`/`LAMBDA_NM` (operating wavelength, 810 nm),
`TAU_GAS_ZENITH`, `AEROSOL_SITES` (per-city AOD/Angstrom parameters).

| Function | Signature | Returns |
|---|---|---|
| `tau_rayleigh_zenith` | `(lam_um, P_site_hpa=1013.25, P_std_hpa=1013.25)` | Rayleigh (molecular) zenith optical depth, pressure-scaled. |
| `tau_aerosol_zenith` | `(aod0, alpha, lam0_nm, lam_nm=LAMBDA_NM)` | Aerosol zenith optical depth via Angstrom power-law wavelength scaling. |
| `kim_kruse_q` | `(V_km)` | Kim/Kruse aerosol size-distribution parameter, selected by visibility regime. |
| `tau_aerosol_visibility` | `(V_km, lam_nm=LAMBDA_NM, H_aerosol_km=2.0)` | Visibility-based aerosol optical depth cross-check; returns `(tau, q)`. |
| `kasten_young_airmass` | `(zenith_deg)` | Kasten & Young (1989) airmass (zenith optical depth -> slant optical depth scale factor). |
| `site_tau_zenith` | `(tau_R_810)` | Dict of total zenith optical depth (Rayleigh + aerosol + gas) per site, baseline case. |
| `site_tau_zenith_winter` | `(tau_R_810)` | Same, Delhi winter-haze sensitivity case. |

---

## `src/atmosphere/turbulence.py` — slant-path turbulence & scintillation

Module constants: `LAMBDA_PHASE4_M`, `K_810`, `H_TOP_M`,
`SPHERICAL_UPLINK_FACTOR`, `D_PROVISIONAL_M`, `TURBULENCE_SITE_CASES`.

| Function | Signature | Returns |
|---|---|---|
| `cn2_hv` | `(h_m, A0=1.7e-14, v=21.0)` | Hufnagel-Valley `Cn2(h)` refractive-index structure parameter. |
| `fried_r0` | `(zenith_deg, lam_m, A0=1.7e-14, v=21.0, h_top=H_TOP_M)` | `(r0, integrated_Cn2)`: Fried parameter and the path-integrated turbulence strength. |
| `sigma_R2_slant_planewave` | `(zenith_deg, lam_m, cn2_profile_fn, h0=0.0, h_top=H_TOP_M)` | General plane-wave (downlink) slant-path Rytov variance for an arbitrary `Cn2(h)` profile. |
| `sigma_R2_slant` | `(zenith_deg, lam_m, cn2_profile_fn, h0=0.0, h_top=H_TOP_M)` | **Corrected spherical-wave (uplink) Rytov variance** — the function used everywhere downstream, since SNS-TF-QKD here is ground -> satellite. |
| `aoa_jitter_urad` | `(zenith_deg, D=D_PROVISIONAL_M, lam_m=LAMBDA_PHASE4_M)` | Angle-of-arrival jitter (Tyler 1994), microradians rms. |
| `strehl_uncompensated` | `(r0, D=D_PROVISIONAL_M)` | Marechal/Noll uncompensated Strehl ratio from `D/r0`. |

---

## `src/atmosphere/seasonal.py` — seasonal datasets

Pure data module: `SEASONAL_AEROSOL`, `A0_BASELINE`, `SEASONAL_TURBULENCE`,
`SEASONS`. No functions; consumed directly by
`src/experiments/seasonal_experiments.py` and downstream stages that need a
season-specific aerosol/turbulence case.

---

## `src/pat/losses.py` — PAT, optics and detector real-loss module

Module constants: hardware baseline (`LAMBDA_P5_M`, `THETA_DIV_RAD`,
`D_RX_M`, `G_TX`, `G_RX`, `ETA_RX_INTERNAL`, `POINTING_CASES`), detector
realism (`ETA_DET_REAL`, `DARK_COUNT_RATE_HZ`, `P_DARK_REAL`), background
light (`R_NOISE_BASELINE_HZ`, `BACKGROUND_SCENARIOS`), `E_A_MISALIGNMENT`.

| Function | Signature | Returns |
|---|---|---|
| `eta_diffraction` | `(L_m, G_tx=G_TX, G_rx=G_RX, lam=LAMBDA_P5_M)` | Friis-equation diffraction/collection efficiency at range `L_m`. |
| `eta_pointing` | `(sigma_jitter_rad, theta_div=THETA_DIV_RAD, theta_bias=0.0)` | Gaussian-beam pointing-loss factor from jitter and (optional) static bias. |
| `optimize_instantaneous_real` | `(eta_value, e_a, p_dark, eps_bounds=(0.01, 0.99), mu_bounds=(0.01, 1.0), n_grid=20)` | Grid + Nelder-Mead `(eps, mu')` optimization with an explicit `p_dark` override; returns `(R_opt, eps_opt, mu_opt)`. |

---

## `src/coupling/phase6.py` — ensemble optimizer with a real p_dark

| Function | Signature | Returns |
|---|---|---|
| `optimize_ensemble_real` | `(eta_det, e_a, sigma_R2, p_dark, n_mc=60_000, seed=0, ...)` | Grid + Nelder-Mead `(eps, mu')` optimization under turbulence-faded (Gamma-Gamma Monte Carlo) `eta`, with an explicit `p_dark`; returns `(R_opt, eps_opt, mu_opt)`. |

---

## `src/weather/availability.py` — climatology-based availability

Module constants: `RAINY_DAYS_PER_MONTH`, `DAYS_IN_MONTH`, `FOG_DAYS_DELHI`.

| Function | Signature | Returns |
|---|---|---|
| `availability_fraction` | `(station, month)` | Fraction of the month expected to be usable (`1 - outage_days/days_in_month`), from rainy-day + (Delhi only) winter-fog climatology. |

---

## `src/optimization/joint.py` — joint signal/decoy optimization

| Function | Signature | Returns |
|---|---|---|
| `sns_key_rate_eta_joint` | `(eps, mu2, mu1, eta, e_a, f_ec=F_EC, p_dark=P_DARK)` | Same physics as `src.engine.sns_key_rate_eta`, but with `mu1` as a free argument instead of a fixed ratio of `mu2`. |
| `refactor_correctness_check` | `()` | Asserts `sns_key_rate_eta_joint` at `mu1=0.2*mu2` reproduces `sns_key_rate_eta` exactly; prints and returns the max deviation. |
| `optimize_joint_signal_decoy` | `(eta_value, e_a, p_dark, n_grid=14, eps_bounds=(1e-5, 0.5), mu2_bounds=(1e-3, 1.0), frac_bounds=(0.005, 0.95))` | Log-grid + Nelder-Mead joint optimization over `(eps, mu2, mu1=frac*mu2)`; returns `(R_opt, eps_opt, mu2_opt, mu1_opt)`. |
| `optimize_joint_ensemble` | `(eta_det, e_a, sigma_R2, p_dark, n_mc=20_000, seed=0, n_grid=10, ...)` | Same, under turbulence-faded (Gamma-Gamma Monte Carlo) `eta`. |

---

## `src/links/` — multi-city-pair links

### `city_pairs.py`

Pure data module: `CITY_DB` (6-city dict: `lat, lon, alt_km, aod500, alpha,
turb_mult, conf`), `LINK_PAIRS` (the 4 evaluated pairs), `LINK_COLORS`.

### `pass_builder.py`

| Function | Signature | Returns |
|---|---|---|
| `build_pass` | `(city_a, city_b, orb, a, inc, gmst0, el_min=20.0, t_half=500.0, n_t=2001)` | Engineers a great-circle-midpoint crossing for any two cities on a given fixed orbital plane; returns a dict of `dist_km, t, el_a, el_b, rng_a, rng_b, common, dt`. |
| `evaluate_link` | `(city_a, city_b, pr, eta_point_val, suppression_factor=1000.0, e_a=E_A_MISALIGNMENT)` | Fully-dressed real-loss channel + `optimize_joint_ensemble` at the link's minimum-joint-slant-range point; returns a dict keyed `"A"`/`"B"` with per-station `eta_total, sigma_R2, R, eps, mu2, mu1`. |

### `metrics.py`

Module constant: `CITY_BG_MULTIPLIER` (per-city night background-light
multiplier).

| Function | Signature | Returns |
|---|---|---|
| `gain_X_asym` | `(mu, eta_A, eta_B, p_dark=P_DARK)` | X-window gain with each arm's transmittance kept separate (additive combination). |
| `error_X_asym` | `(mu, eta_A, eta_B, e_a, p_dark=P_DARK)` | X-window QBER, asymmetric-arm version. |
| `S_Z_gain_asym` | `(eps, mu2, eta_A, eta_B, p_dark=P_DARK)` | Z-window gain, split into the four true cases (Alice-only / Bob-only / both / neither send). |
| `E_Z_qber_asym` | `(eps, mu2, eta_A, eta_B, p_dark=P_DARK)` | Z-window QBER, asymmetric-arm version. |
| `sns_key_rate_asym_full` | `(eps, mu2, mu1, eta_A, eta_B, e_a, f_ec=F_EC, p_dark=P_DARK)` | Full asymmetric rate-formula dict (`S1,S2,E1,s0,s1,e1ph,Sz,Ez,R`). |
| `validate_symmetric_limit` | `(eta_t=0.01, p_dark=P_DARK)` | Asserts the asymmetric formulas reduce exactly to the symmetric `src.engine` formulas at `eta_A=eta_B`; returns the max deviation. |
| `optimize_joint_asym_ensemble` | `(eta_A_det, eta_B_det, sigma_A, sigma_B, e_a, p_dark, n_mc=4000, n_grid=6, seed=0, ...)` | Joint `(eps, mu2, mu1)` optimization with **independent** per-arm Gamma-Gamma fading; returns `(R_opt, eps_opt, mu2_opt, mu1_opt)`. |

---

## `src/experiments/` — stage orchestrators

Each file exposes one `run_full_partN_pipeline(...)` entry point per stage.
These functions run that stage's computation in order (printing tables and
producing figures as they go) and return a dict of that stage's key outputs,
which the next stage's orchestrator takes as an explicit argument. This is
the layer `scripts/run_part3_to_13_satellite.py` calls into; see
`docs/workflow.md` for the full call graph and what each stage returns.

| File | Entry point |
|---|---|
| `satellite_experiments.py` | `run_full_part3_pipeline(show_plots=True)` |
| `atmosphere_experiments.py` | `run_full_part4_pipeline(state, show_plots=True)` |
| `turbulence_experiments.py` | `run_full_part5_pipeline(state, show_plots=True)` |
| `seasonal_experiments.py` | `run_full_part6_pipeline(state, tau_R_810, show_plots=True)` |
| `pat_experiments.py` | `run_full_part7_pipeline(state, atm_results, show_plots=True)` |
| `phase6_experiments.py` | `run_full_part8_pipeline(state, atm_results, turb_results_site, phase5_results, show_plots=True)` |
| `weather_experiments.py` | `run_full_part9_pipeline(state, atm_state, turb_state, phase5_state, phase6_state, show_plots=True)` |
| `phase10_experiments.py` | `run_full_part10_pipeline(state, phase5_state, turb_state, phase6_state, weather_state, show_plots=True)` |
| `consolidation_experiments.py` | `run_full_part11_pipeline(state, atm_state, turb_state, phase5_state, phase6_state, weather_state, phase10_state, show_plots=True)` |
| `city_links_experiments.py` | `run_full_part12_pipeline(state, show_plots=True)` |
| `link_metrics_experiments.py` | `run_full_part13_pipeline(links_state, show_plots=True)` |
| `fiber_experiments.py` | fiber-channel demos, tables, and figures |
| `fso_experiments.py` | `run_full_part2_pipeline(...)` — full FSO pipeline |
