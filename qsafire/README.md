# SNS-TF-QKD Simulation Suite

*Developed as part of a Summer Internship at the Scientific Analysis Group
(SAG), DRDO.*

A modular Python package for simulating and analyzing **Sending-or-Not-Sending
Twin-Field Quantum Key Distribution (SNS-TF-QKD)** across three physical
channel types:

1. **Fiber** — a standard telecom-fiber point-to-point link.
2. **Free-space-optical (FSO) with atmospheric turbulence** — a horizontal
   or near-ground FSO link with Gamma-Gamma scintillation.
3. **Ground-to-satellite optical links** — a full physical-layer simulation
   of an uplink to a sun-synchronous satellite, including orbital mechanics,
   atmospheric transmittance, slant-path turbulence, seasonal variation,
   realistic pointing/acquisition/tracking (PAT) and detector hardware
   losses, weather-gated availability, joint intensity/probability
   optimization, and multi-city-pair link comparison.

The package implements the asymptotic key-rate formula of:

> X.-B. Wang, Z.-W. Yu, X.-L. Hu, *"Sending or not sending: Twin-field
> quantum key distribution with large misalignment error"*, Phys. Rev. A
> **98**, 062323 (2018), arXiv:1805.09222,

and extends it with a from-first-principles free-space channel model,
realistic hardware/environment loss budgets, and a fully asymmetric
(per-arm) generalization of the protocol's rate formulas for links where the
two ground stations do not share the same channel transmittance.

## Features

- **Exact protocol fidelity.** The core SNS-TF-QKD equations (photon-number
  weights, gain/error model, Z-window gain/QBER, decoy-state yield and
  phase-error bounds, the asymptotic key-rate formula) live in one canonical
  module (`src/engine.py`) and are reused, unmodified, by every channel
  model in the package — fiber, FSO, and satellite.
- **Physically-grounded free-space channel model.** Rytov-variance
  turbulence strength, Gamma-Gamma irradiance fading, and turbulence-regime
  classification, cross-validated between a closed-form Gauss-Laguerre
  quadrature evaluation and direct Monte Carlo sampling.
- **A complete ground-to-satellite link budget**, built up in physically
  separable layers:
  - WGS-84 geodesy and a *solved* (not assumed) sun-synchronous orbit design
  - Rayleigh + aerosol + gaseous zenith optical depth, scaled to the slant
    path via the Kasten & Young airmass formula
  - Hufnagel-Valley Cn2(h) turbulence profile, Fried parameter, and the
    corrected spherical-wave (uplink) Rytov variance
  - an explicit winter/summer/monsoon seasonal extension
  - real transmitter/receiver hardware: Friis-equation diffraction loss,
    Gaussian pointing-jitter loss, a characterized avalanche photodiode,
    and day/night background-light scenarios
  - Monte Carlo ensemble optimization that folds turbulence fading into the
    real-loss channel and derives the background-suppression margin a real
    system would need
  - weather-gated annual availability from real rainy-day/fog climatology
  - joint (not fixed-ratio) signal/decoy intensity and sending-probability
    optimization
  - a six-city link database and a fully asymmetric per-arm generalization
    of the protocol's rate formulas for links between two different cities
- **Every stage is independently validated** against an external reference
  (a published formula, a textbook target, a real reported link budget, or
  an internal algebraic identity) rather than only checked for
  self-consistency — see `docs/theory.md` for the full list.
- **No duplicated physics.** Every module that needs the SNS-TF-QKD rate
  equations imports them from `src/engine.py`; every module that needs
  Gamma-Gamma fading imports it from `src/fso/channel_models.py`. Downstream
  layers (satellite, links) build on top of the validated fiber/FSO core
  instead of re-deriving it.
- **Figure-quality plots** for every stage: pass geometry, atmospheric and
  turbulence time series, key rate vs. time/elevation/season, sensitivity
  tornado charts, and multi-link comparison figures.
- **CSV / figure / text-report export** for the fiber+FSO pipeline.

## Project Structure

```
project/
├── README.md
├── requirements.txt
├── pyproject.toml
├── config/
│   └── config.py                  # output-directory paths
├── src/
│   ├── engine.py                  # canonical SNS-TF-QKD rate equations (eta-native)
│   ├── simulation/
│   │   └── fiber_channel.py       # fiber channel model + optimizer
│   ├── fso/
│   │   ├── config.py               # FSO shared constants
│   │   ├── channel_models.py       # deterministic FSO loss + Gamma-Gamma turbulence
│   │   ├── monte_carlo.py          # Monte Carlo ensemble averaging
│   │   ├── quadrature.py           # Gauss-Laguerre quadrature ensemble averaging
│   │   ├── optimization.py         # instantaneous / ensemble / joint optimization
│   │   ├── metrics.py              # performance-metrics table builders
│   │   ├── visualization.py        # publication figures
│   │   └── export.py               # report + file export
│   ├── satellite/
│   │   └── geometry.py             # WGS-84 geodesy, orbit design, propagation, look angles
│   ├── atmosphere/
│   │   ├── transmittance.py        # Rayleigh + aerosol + gas zenith optical depth
│   │   ├── turbulence.py           # Hufnagel-Valley profile, Fried r0, slant Rytov variance
│   │   └── seasonal.py             # winter/summer/monsoon aerosol + turbulence datasets
│   ├── pat/
│   │   └── losses.py               # diffraction, pointing, detector, background-light losses
│   ├── coupling/
│   │   └── phase6.py               # ensemble optimizer with a real (suppressed) p_dark
│   ├── weather/
│   │   └── availability.py         # rainy-day / fog climatology -> availability fraction
│   ├── optimization/
│   │   └── joint.py                # joint signal+decoy intensity optimizer
│   ├── links/
│   │   ├── city_pairs.py           # 6-city parameter database
│   │   ├── pass_builder.py         # generalized pass construction + link evaluation
│   │   └── metrics.py              # fully asymmetric (per-arm) rate formulas
│   ├── visualization/
│   │   └── fiber_plots.py          # fiber-channel plots
│   └── experiments/                # orchestration layer: one run_full_pipeline() per stage,
│       └── *.py                    #   printing tables and producing figures in sequence
├── scripts/
│   ├── train.py                        # fiber + FSO pipeline
│   ├── run_part1_fiber.py              # fiber channel only
│   ├── run_part2_fso.py                # FSO/turbulence channel only
│   └── run_part3_to_13_satellite.py    # full satellite-link pipeline
├── dashboard/
│   ├── app.py                      # Streamlit dashboard over the full satellite pipeline
│   └── README.md
├── tests/
│   └── smoke_test.py               # numerical-regression checks
└── docs/
    ├── architecture.md             # module dependency graph
    ├── workflow.md                 # execution-flow walkthrough
    ├── module_documentation.md     # per-module purpose/inputs/outputs
    ├── api.md                      # full function-level API reference
    └── theory.md                   # physics/math background, module by module
```

## Installation

```bash
python -m venv .venv
source .venv/bin/activate          # on Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Requirements

- Python >= 3.9
- numpy, pandas, matplotlib, scipy, tqdm (see `requirements.txt`)

## Configuration

Physics constants (detector efficiency, dark-count probability,
error-correction factor, decoy ratio, fiber/FSO attenuation, turbulence Cn2
values, hardware parameters, distance grids, optimization bounds) live
directly next to the equations that use them, in the module where they are
physically relevant (`src/engine.py`, `src/fso/config.py`,
`src/pat/losses.py`, `src/atmosphere/*.py`), so that reorganizing modules
never risks silently changing a numeric value used in a computation.
`config/config.py` only centralizes output-directory paths.

## Usage

### Running an entire pipeline

```bash
python scripts/train.py                       # fiber + FSO, with plots
python scripts/train.py --no-plots            # same, no matplotlib rendering (faster, CI-friendly)
python scripts/run_part3_to_13_satellite.py            # full satellite-link pipeline, with plots
python scripts/run_part3_to_13_satellite.py --no-plots # same, no figures (faster)
```

The satellite pipeline runs Monte Carlo ensemble optimization at many points
along the pass and across four city-pair links, so
`run_part3_to_13_satellite.py --no-plots` takes on the order of a few
minutes end to end — the geometry/atmosphere/turbulence stages are each
fast; the cumulative Monte Carlo cost in the later optimization stages
dominates the runtime.

### Interactive dashboard

```bash
streamlit run dashboard/app.py
```

Opens a browser dashboard over the full satellite-link pipeline — pass
geometry, atmosphere, turbulence, real-loss key rate, weather-gated annual
availability, multi-city links, and the final per-link metrics table, all in
one place. First load takes a few minutes (cached afterward); see
`dashboard/README.md`.

### Running a single stage

```bash
python scripts/run_part1_fiber.py    # fiber SNS-TF-QKD only
python scripts/run_part2_fso.py      # FSO/turbulence channel only
```

### Using the library directly

```python
from src.simulation.fiber_channel import sns_key_rate, optimize_key_rate

R = sns_key_rate(eps=0.05, mu_prime=0.3, L=200, e_a=0.15)
R_opt, eps_opt, mu_opt = optimize_key_rate(L=300, e_a=0.25)
```

```python
from src.fso.quadrature import gauss_laguerre_ensemble_rate
from src.fso.channel_models import rytov_variance
from src.fso.config import CN2_MODERATE

sigma2 = rytov_variance(d_km=6.0, Cn2=CN2_MODERATE)
R_ensemble = gauss_laguerre_ensemble_rate(eps=0.02, mu_prime=0.3, d_km=6.0,
                                           e_a=0.15, sigma_R2=sigma2)
```

```python
from src.experiments.satellite_experiments import run_full_part3_pipeline
from src.experiments.atmosphere_experiments import run_full_part4_pipeline

state = run_full_part3_pipeline(show_plots=False)       # pass geometry
atm = run_full_part4_pipeline(state, show_plots=False)  # atmospheric transmittance
```

See `docs/api.md` for the full function reference and `docs/theory.md` for
the physics behind each stage.

## Outputs

Running the fiber/FSO scripts writes:

- `outputs/csv/*.csv` — deterministic FSO table, Gamma-Gamma table,
  instantaneous-optimization map, maximum-performance summary.
- `outputs/csv/summary_report.txt` — plain-text summary report.
- `outputs/figures/*.png` and `*.pdf` — publication figures.

The satellite pipeline writes its figures (`phase*.png`) to the working
directory and prints its tables directly; see `docs/workflow.md` for what
each stage produces.

## Troubleshooting

- **The FSO Gamma-Gamma table is slow**: `build_gamma_gamma_table` runs
  Differential Evolution + Basin Hopping + grid search at every distance and
  turbulence case (16 distances x 3 cases by default) — this is
  compute-heavy by design (triangulated global optimization). Use
  `scripts/run_part2_fso.py --no-plots` and/or reduce
  `src/fso/config.DISTANCES_MAIN` for faster iteration during development.
- **The satellite pipeline is slow**: the coupling, joint-optimization, and
  per-link metrics stages (`src/coupling`, `src/optimization`,
  `src/links/metrics.py`) run Monte Carlo ensemble optimization; reduce
  `n_mc`/`n_grid` arguments where exposed for faster iteration.
- **`ModuleNotFoundError: No module named 'src'`**: run scripts from the
  project root, or `pip install -e .` — the provided `pyproject.toml`
  supports `pip install -e .` directly.
- **Headless / CI environments**: set `MPLBACKEND=Agg` before running, or
  pass `--no-plots`.

## Future Improvements

- Replace the plain-text smoke tests with a `pytest` suite.
- Re-derive the background-suppression factor per city-pair link instead of
  reusing the Delhi-Mumbai-derived value uniformly (`src/links/pass_builder.py`).
- Extend the seasonal/weather-gated annual analysis to the three additional
  city-pair links (currently only Delhi-Mumbai has a full annual model).
- Vectorize the per-distance optimization loop in
  `build_gamma_gamma_table` for faster full-grid runs.

## Acknowledgment

This project was carried out as part of a Summer Internship at the
Scientific Analysis Group (SAG), Defence Research and Development
Organisation (DRDO).

