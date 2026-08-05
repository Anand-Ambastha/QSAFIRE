# Module Documentation

A quick-reference index of every module's purpose, inputs/outputs, and
dependencies. For full function signatures see `docs/api.md`; for the
physics/math behind each module see `docs/theory.md`.

## `config/config.py`
**Purpose:** centralizes output-directory paths (`outputs/figures`,
`outputs/csv`). **Inputs:** none. **Outputs:** `PROJECT_ROOT`,
`OUTPUTS_DIR`, `FIGURES_DIR`, `CSV_DIR`. **Dependencies:** `os`.

## `src/engine.py`
**Purpose:** canonical, eta-native SNS-TF-QKD physics engine — photon-number
weights, gain/error model, Z-window gain/QBER, decoy-state single-photon
yield and phase-error bounds, binary entropy, the asymptotic key-rate
formula. Every other channel model in the package imports from here rather
than re-deriving the protocol. **Inputs:** `mu`, `eta`, `eps`, `e_a`, and
protocol constants (`P_DARK`, `F_EC`, `DECOY_RATIO`). **Outputs:** scalars
or numpy arrays (fully vectorized over `eta`). **Dependencies:** `numpy`.

## `src/simulation/fiber_channel.py`
**Purpose:** fiber channel model (`channel_eta(L)`), the distance-driven
key-rate wrapper (`sns_key_rate`), and its two-stage (grid + local)
optimizer. **Dependencies:** `numpy`, `scipy.optimize`, `src.engine`.

## `src/fso/config.py`
**Purpose:** shared FSO-package constants — operating wavelength,
attenuation coefficient, distance grids, turbulence-case definitions
(`CN2_WEAK/MODERATE/STRONG`), optimization bounds, Monte Carlo defaults,
output paths. No functions.

## `src/fso/channel_models.py`
**Purpose:** deterministic FSO loss, the horizontal-path Rytov variance,
Gamma-Gamma turbulence-fading parameters and sampling, turbulence-regime
classification. **Dependencies:** `numpy`, `scipy.special`.

## `src/fso/monte_carlo.py`
**Purpose:** Monte Carlo ensemble averaging of the key rate under
Gamma-Gamma fading, plus a convergence-study helper. **Dependencies:**
`numpy`, `src.fso.channel_models`, `src.engine`.

## `src/fso/quadrature.py`
**Purpose:** Gauss-Laguerre quadrature ensemble averaging — a closed-form
alternative to Monte Carlo for repeated evaluation inside an optimizer, with
a validation helper that cross-checks both against each other.
**Dependencies:** `numpy`, `scipy.special`, `src.fso.channel_models`,
`src.engine`.

## `src/fso/optimization.py`
**Purpose:** instantaneous (deterministic-channel) and ensemble
(turbulence-faded) optimization of `(eps, mu')`, triangulating local grid
search with global Differential Evolution and Basin Hopping.
**Dependencies:** `numpy`, `scipy.optimize`, `src.fso.quadrature`,
`src.fso.channel_models`.

## `src/fso/metrics.py`
**Purpose:** builds the deterministic-channel and Gamma-Gamma performance
tables (DataFrames) across a distance grid and turbulence cases, and the
maximum-performance summary. **Dependencies:** `pandas`,
`src.fso.optimization`, `src.fso.channel_models`.

## `src/fso/visualization.py`, `src/fso/export.py`
**Purpose:** publication figures and CSV/text-report export for the FSO
pipeline, called from `src/experiments/fso_experiments.py`.
**Dependencies:** `matplotlib`, `pandas`.

## `src/satellite/geometry.py`
**Purpose:** WGS-84 geodetic<->ECEF conversion, sun-synchronous orbit design
(solved via the J2 secular nodal-precession condition, not assumed), GMST
Earth orientation, ECI/ECEF propagation, ENU topocentric look angles
(elevation/azimuth/slant-range), great-circle geometry. **Dependencies:**
`numpy` only.

## `src/atmosphere/transmittance.py`
**Purpose:** Beer's-law static zenith optical depth (Rayleigh + aerosol +
gas) and Kasten & Young (1989) airmass scaling to the slant path.
**Dependencies:** `numpy`.

## `src/atmosphere/turbulence.py`
**Purpose:** Hufnagel-Valley Cn2(h) profile, Fried parameter, and the
corrected spherical-wave (uplink) Rytov variance as path integrals;
provisional angle-of-arrival jitter and Strehl ratio. **Dependencies:**
`numpy`, `scipy.integrate`.

## `src/atmosphere/seasonal.py`
**Purpose:** cited winter/summer/monsoon aerosol and ground-turbulence
datasets, per city. Pure data module, no functions.

## `src/pat/losses.py`
**Purpose:** real hardware losses — Friis-equation diffraction/collection
efficiency, Gaussian-beam pointing loss, real Si-APD detector parameters,
day/night + per-city background-light scenarios, an instantaneous optimizer
that threads a real `p_dark` through `src.engine.sns_key_rate_eta`.
**Dependencies:** `numpy`, `scipy.optimize`, `src.engine`.

## `src/coupling/phase6.py`
**Purpose:** Monte Carlo ensemble (turbulence-faded) optimizer with a real
`p_dark` pass-through, reusing `gamma_gamma_params`/`sample_gamma_gamma`
unchanged from `src.fso.channel_models`. **Dependencies:** `numpy`,
`scipy.optimize`, `src.engine`, `src.fso.channel_models`.

## `src/weather/availability.py`
**Purpose:** real rainy-day-count climatology (plus a winter-fog proxy for
Delhi), combined into a monthly cloud-free-availability fraction. Pure data
+ one function, no external dependencies beyond the module itself.

## `src/optimization/joint.py`
**Purpose:** corrects the epsilon search bound to the protocol-correct `(0,
0.5)` and independently optimizes the decoy intensity instead of a fixed
ratio; includes a refactor-correctness check against `src.engine`.
**Dependencies:** `numpy`, `scipy.optimize`, `src.engine`,
`src.fso.channel_models`.

## `src/links/city_pairs.py`
**Purpose:** a 6-city parameter database (location, altitude, aerosol,
turbulence multiplier) and the 4 evaluated link pairs. Pure data module.

## `src/links/pass_builder.py`
**Purpose:** generalized pass construction (any two cities, on a fixed
orbital plane) and per-link real-loss channel evaluation, reusing every
`src.satellite`/`src.atmosphere`/`src.pat`/`src.optimization` function
unchanged. **Dependencies:** `numpy`, `src.satellite.geometry`,
`src.atmosphere.*`, `src.pat.losses`, `src.optimization.joint`,
`src.links.city_pairs`.

## `src/links/metrics.py`
**Purpose:** the fully asymmetric (per-arm) SNS-TF-QKD rate formulas —
validated to reduce exactly to `src.engine`'s symmetric formulas when
`eta_A = eta_B` — plus the independent-per-arm ensemble optimizer.
**Dependencies:** `numpy`, `scipy.optimize`, `src.engine`,
`src.fso.channel_models`.

## `src/experiments/*`
**Purpose:** orchestration layer — one `run_full_pipeline(...)` per stage,
running that stage's computation in order (tables, figures) and returning a
dict of key outputs for the next stage. See `docs/api.md` for each stage's
entry point and `docs/workflow.md` for the full call graph.

## `src/visualization/fiber_plots.py`
**Purpose:** fiber-channel plots (key rate vs. distance, sensitivity
sweeps). **Dependencies:** `matplotlib`.
