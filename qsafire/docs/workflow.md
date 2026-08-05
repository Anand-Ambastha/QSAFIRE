# Execution Workflow

## Fiber channel (`scripts/run_part1_fiber.py`)

Driven by `src/experiments/fiber_experiments.py`:

1. `demo_channel_scaling()` — prints the `eta(L)` distance-scaling sanity
   table.
2. `demo_gain_and_error()` — plots gain/error vs. decoy intensity at a
   representative distance.
3. `demo_z_qber()` — plots Z-basis QBER vs. the sending probability `eps`.
4. `validate_s1_decoy_estimate()` — compares the decoy-state lower bound on
   `s1` against the true single-photon yield across distances.
5. `demo_e1ph_recovery()` — recovers the phase-error upper bound vs.
   injected misalignment error.
6. `quick_check_key_rate()` / `quick_check_optimizer()` — single-point
   sanity checks of the key-rate formula and the `(eps, mu')` optimizer.
7. `figure1_sweep()` — key rate vs. distance at several misalignment levels,
   each point independently optimized.
8. `figure2_sweep()` — key rate vs. misalignment error at a fixed distance.
9. `validation_tables()` — full diagnostic table of every intermediate
   quantity at representative (distance, misalignment) pairs.

## FSO + Gamma-Gamma turbulence (`scripts/run_part2_fso.py`)

Driven by `src/experiments/fso_experiments.run_full_part2_pipeline`:

1. **Channel and regime validation** — validates the Gamma-Gamma sampler
   against theoretical variance/PDF normalization, then plots Rytov
   variance and the turbulence-regime map vs. distance.
2. **Monte Carlo / quadrature cross-validation** — cross-validates the
   Gauss-Laguerre quadrature ensemble-rate estimator against a 150k-sample
   Monte Carlo run at every distance (target: <1% relative error), then
   runs a Monte Carlo convergence diagnostic.
3. **Full performance pipeline**:
   - Builds the deterministic-FSO and Gamma-Gamma performance-metrics
     tables via `optimize_instantaneous` and `joint_optimize_ensemble`
     (grid search + Differential Evolution + Basin Hopping, keeping the
     best result found).
   - Generates all publication figures: SKR/QBER/error/yield/loss vs.
     distance for both channel models, optimal-parameter curves, the
     epsilon-improvement figure, instantaneous alpha/epsilon maps, and the
     joint-optimization landscape heatmap.
   - Computes the maximum-performance summary (max SKR, max secure
     distance, max tolerable loss per channel model).
4. **Export** — writes the plain-text summary report and lists every
   exported CSV/figure file.

## Fiber + FSO combined (`scripts/train.py`)

Runs the fiber pipeline in full, then the FSO pipeline in full. Pass
`--no-plots` to skip matplotlib rendering (useful for headless/CI runs);
all CSV/report outputs are still written.

```
config (constants) -> engine.py / fso/*  -> experiments/* (orchestration,
                                              prints, dataframes)
                                          -> visualization (figures, PNG/PDF)
                                          -> fso/export.py (CSV, .txt report)
```

All state is passed explicitly through function arguments and return
values (dataframes, dicts, numpy arrays) — no hidden global state.

## Satellite link budget (`scripts/run_part3_to_13_satellite.py`)

Runs the full ground-to-satellite pipeline end to end. Each stage's
`run_full_partN_pipeline(...)` function takes the previous stages' returned
state dict(s) as explicit arguments and returns its own dict, so the full
chain is:

```
run_full_part3_pipeline()              -> state          (pass geometry)
run_full_part4_pipeline(state)         -> atm_state       (atmospheric transmittance)
run_full_part5_pipeline(state)         -> turb_state       (slant-path turbulence)
run_full_part6_pipeline(state, ...)    -> seasonal_state   (seasonal extension)
run_full_part7_pipeline(state, ...)    -> phase5_state      (real hardware losses)
run_full_part8_pipeline(state, ...)    -> phase6_state       (ensemble coupling)
run_full_part9_pipeline(state, ...)    -> weather_state       (weather-gated availability)
run_full_part10_pipeline(state, ...)   -> phase10_state        (joint signal/decoy optimization)
run_full_part11_pipeline(state, ...)   -> phase11_state         (validation consolidation)
run_full_part12_pipeline(state)        -> links_state             (multi-city links)
run_full_part13_pipeline(links_state)  -> link_metrics_state       (asymmetric per-link metrics)
```

Every dependency between stages is an explicit function argument rather
than shared global state, so any stage can be re-run standalone given the
prior stage's returned dict (useful for iterating on one physical effect
without re-running the whole chain — see `docs/api.md` for each stage's
exact entry-point signature).

The geometry, atmosphere, and turbulence stages each run in well under a
second; the coupling, joint-optimization, and per-link-metrics stages run
Monte Carlo ensemble optimization at many points along the pass and across
four city-pair links, and dominate the pipeline's total runtime (a few
minutes end to end with `--no-plots`).
