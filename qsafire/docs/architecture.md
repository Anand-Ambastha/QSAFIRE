# Architecture

## Overview

The package is organized as a set of physically-separable layers, all
built on top of one shared protocol core. Each layer only depends on the
layers below it, and every layer's numeric output is independently
validated (see `docs/theory.md`) rather than only checked for internal
self-consistency.

```
                       ┌───────────────────────┐
                       │    src/engine.py       │  canonical eta-native
                       │  (symmetric SNS-TF-QKD)│  SNS-TF-QKD physics
                       └───────────┬────────────┘
                                   │
        ┌──────────────┬──────────┴──────────┬───────────────────┐
        │               │                      │                   │
┌───────▼────────┐ ┌────▼─────────────┐ ┌──────▼───────────┐ ┌─────▼──────────┐
│ src/simulation/  │ │ src/fso/          │ │ src/pat/          │ │ src/links/      │
│  fiber_channel.py│ │  channel_models,  │ │  losses.py        │ │  metrics.py      │
│  (fiber eta(L))  │ │  monte_carlo,     │ │ (real hardware    │ │ (asymmetric per- │
│                  │ │  quadrature,      │ │  loss budget)      │ │  arm generaliz-  │
│                  │ │  optimization,    │ │                    │ │  ation, reduces  │
│                  │ │  metrics          │ │                    │ │  exactly to      │
│                  │ │ (Gamma-Gamma      │ │                    │ │  engine.py at    │
│                  │ │  turbulence)      │ │                    │ │  eta_A=eta_B)    │
└──────────────────┘ └───────┬───────────┘ └──────────┬─────────┘ └────────┬─────────┘
                              │                        │                    │
              ┌───────────────┴────────────────────────┴────────────────────┘
              │
   ┌──────────▼───────────┐   ┌───────────────────┐   ┌────────────────────┐
   │ src/satellite/         │   │ src/atmosphere/     │   │ src/weather/         │
   │  geometry.py           │──▶│  transmittance.py   │   │  availability.py     │
   │ (orbit, pass geometry) │   │  turbulence.py      │   │ (rainy-day/fog       │
   │                        │   │  seasonal.py        │   │  climatology)        │
   └────────────┬───────────┘   └──────────┬──────────┘   └──────────┬───────────┘
                │                          │                          │
                └──────────────┬───────────┴──────────────┬───────────┘
                               │                            │
                    ┌──────────▼──────────┐      ┌──────────▼───────────┐
                    │ src/coupling/        │      │ src/optimization/     │
                    │  phase6.py            │      │  joint.py              │
                    │ (ensemble optimizer,  │      │ (joint signal+decoy   │
                    │  real p_dark)         │      │  optimization)         │
                    └──────────┬────────────┘      └──────────┬────────────┘
                               │                               │
                               └───────────────┬───────────────┘
                                               │
                                  ┌─────────────▼──────────────┐
                                  │  src/links/pass_builder.py  │
                                  │  (any two cities, reuses    │
                                  │   every layer above)        │
                                  └─────────────────────────────┘

                  All layers are orchestrated by src/experiments/*.py
                  and driven by scripts/*.py.
```

## Design principle: one canonical protocol core

`src/engine.py` implements the SNS-TF-QKD protocol's photon-number weights,
gain/error model, decoy-state bounds, and asymptotic key-rate formula
**once**, in eta-native form (every function takes channel transmittance
`eta` directly, rather than a channel-specific parameter like distance).
Every channel model in the package — fiber (`channel_eta(L)` wraps it),
FSO (deterministic or Gamma-Gamma-faded `eta` drives it), and satellite
(the fully-dressed real-loss `eta` drives it) — calls into this one module.
`src/optimization/joint.py` and `src/links/metrics.py` extend it (decoy
intensity as a free parameter; a fully asymmetric two-arm generalization,
respectively) rather than duplicating it, and both are validated to reduce
exactly to `src.engine`'s output in the relevant special case.

## Design principle: physically separable layers

The satellite link budget is built as a stack of independent physical
effects — geometry, atmosphere, turbulence, season, hardware, background
light, weather — each in its own module, each independently validated
against an external reference where one exists (see `docs/theory.md`).
This means:

- A change to one physical effect (e.g. a better aerosol model) only
  touches one module.
- Every stage's output can be sanity-checked in isolation before it feeds
  the next stage.
- The multi-city-link extension (`src/links/`) can reuse every lower layer
  unchanged — it only adds a city database, generalized pass construction,
  and the asymmetric rate-formula generalization at the top.

## Orchestration layer

`src/experiments/*.py` holds one `run_full_pipeline(...)` function per
stage. Each function runs that stage's computation in order (prints tables,
saves figures) and returns a dict of the stage's key outputs; downstream
stages take the relevant prior dict(s) as explicit function arguments. This
makes every dependency between stages visible in the function signature
rather than implicit in shared global state — see `docs/workflow.md` for
the full call graph. `scripts/*.py` are thin CLI entry points around this
orchestration layer.

## Module responsibilities

| Module | Responsibility |
|---|---|
| `src/engine.py` | Core SNS-TF-QKD equations, eta-native, symmetric two-arm case |
| `src/simulation/fiber_channel.py` | Fiber channel model, distance-driven key rate, optimizer |
| `src/fso/channel_models.py` | Deterministic FSO loss + Gamma-Gamma turbulence-fading model |
| `src/fso/monte_carlo.py`, `quadrature.py` | Two independent ensemble-averaging methods, cross-validated |
| `src/fso/optimization.py`, `metrics.py` | Global optimization + performance tables for the FSO channel |
| `src/satellite/geometry.py` | Geodesy, orbit design, propagation, pass geometry |
| `src/atmosphere/transmittance.py` | Zenith optical depth -> slant-path atmospheric transmittance |
| `src/atmosphere/turbulence.py` | Altitude-resolved turbulence -> slant-path Rytov variance |
| `src/atmosphere/seasonal.py` | Winter/summer/monsoon parameter datasets |
| `src/pat/losses.py` | Real transmitter/receiver/detector/background-light loss budget |
| `src/coupling/phase6.py` | Ensemble optimizer coupling turbulence fading + real losses |
| `src/weather/availability.py` | Climatology-based cloud-free availability |
| `src/optimization/joint.py` | Joint signal/decoy-intensity optimization |
| `src/links/*` | Multi-city database, generalized pass construction, asymmetric protocol |
| `src/experiments/*` | Per-stage orchestration |
