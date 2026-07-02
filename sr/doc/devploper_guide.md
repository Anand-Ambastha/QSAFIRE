# snstfqkd — Developer Reference

**Audience:** engineers integrating with, extending, or maintaining this package.
**Assumes:** familiarity with Python, numpy, and (loosely) the SNS-TF-QKD protocol (Wang, Yu, Hu 2018, arXiv:1805.09222v9). Physics notation is explained inline where it affects API shape, not derived from scratch — see `ARCHITECTURE.md` for the full protocol background and refactor rationale.

---

## Table of Contents

1. [Install & quickstart](#1-install--quickstart)
2. [Package layout](#2-package-layout)
3. [Core API: `channel/`](#3-core-api-channel)
4. [Core API: `protocol/`](#4-core-api-protocol)
5. [Core API: `security/`](#5-core-api-security)
6. [Core API: `optimization/`](#6-core-api-optimization)
7. [`config/` — parameters & presets](#7-config--parameters--presets)
8. [`visualization/` — figure functions](#8-visualization--figure-functions)
9. [`reports/` — CSV & summary reports](#9-reports--csv--summary-reports)
10. [`workflows/` and `cli.py`](#10-workflows-and-clipy)
11. [Extension recipes](#11-extension-recipes)
12. [Testing conventions](#12-testing-conventions)
13. [Performance notes & gotchas](#13-performance-notes--gotchas)
14. [Coding conventions for contributors](#14-coding-conventions-for-contributors)

---

## 1. Install & Quickstart

```bash
pip install -e .
# or: uv pip install -e .
```

Minimal working example — compute a key rate at 300 km:

```python
from snstfqkd.channel import AnalyticalChannel, ChannelParameters
from snstfqkd.security import DecoyStateAnalyzer, SNSKeyRateCalculator

params = ChannelParameters(misalignment_error=0.15)
channel = AnalyticalChannel(params)

analyzer = DecoyStateAnalyzer()
calc = SNSKeyRateCalculator(f=1.1)

mu1, mu2, mu_signal, epsilon, distance = 0.01, 0.15, 0.05, 0.05, 300.0

decoy = analyzer.three_intensity_decoy(
    mu1, mu2,
    channel.vacuum_yield(),
    channel.intensity_yield(mu1, distance),
    channel.intensity_yield(mu2, distance),
    channel.qber_x_basis(mu1, distance),
)
rate = calc.key_rate(
    epsilon, mu_signal, decoy.s1, decoy.e1_ph,
    channel.z_window_yield(mu_signal, distance),
    channel.qber_z_basis(),
)
print(rate)   # 1.3739152435016688e-06
```

Or run the CLI directly:

```bash
uv run -m snstfqkd.main --mode analytical
```

---

## 2. Package Layout

```
snstfqkd/
├── config/          defaults.py, parameters.py, presets.py
├── channel/         base.py, analytical_channel.py, mc_channel.py, detector.py
├── protocol/        statistics.py, phase_slice.py, events.py, pulse.py, protocol_simulator.py
├── security/         decoy_analysis.py, key_rate.py, bounds.py
├── optimization/      optimizer.py
├── visualization/
│   ├── analytical/   key_rate_figs.py, decoy_figs.py, channel_figs.py, optimization_figs.py
│   ├── mc/            core_figs.py, diagnostic_figs.py, performance_figs.py
│   ├── protocol/       window_figs.py, convergence_figs.py
│   └── comparison/      comparison_figs.py
├── workflows/         analytical_workflow.py, mc_workflow.py, protocol_workflow.py, comparison_workflow.py
├── reports/            report_generator.py
├── cli.py, main.py
tests/                  test_*.py (pytest)
```

Every subpackage has an `__init__.py` that re-exports its public symbols — always prefer importing from the subpackage (`from snstfqkd.channel import AnalyticalChannel`) rather than the submodule, unless you need something private (leading underscore).

---

## 3. Core API: `channel/`

### `BaseChannel` (abstract)

All channel implementations satisfy this interface (`snstfqkd/channel/base.py`):

```python
class BaseChannel(ABC):
    params: ChannelParameters

    def channel_transmittance(self, distance: float) -> float: ...
    def total_transmittance(self, distance: float) -> float: ...
    def vacuum_yield(self) -> float: ...
    def intensity_yield(self, intensity: float, distance: float) -> float: ...
    def z_window_yield(self, intensity: float, distance: float) -> float: ...
    def gain(self, intensity: float, distance: float) -> float: ...
    def qber_x_basis(self, intensity: float, distance: float) -> float: ...
    def qber_z_basis(self) -> float: ...
```

Write functions against `BaseChannel`, not a concrete class, so they work for `AnalyticalChannel` and `MonteCarloChannel` interchangeably. This is exactly how `security/`, `optimization/`, and `visualization/analytical`+`mc` are written.

### `ChannelParameters`

```python
@dataclass
class ChannelParameters:
    fiber_loss: float = 0.2                    # dB/km
    detector_efficiency: float = 0.8           # eta_d, (0, 1]
    dark_count_rate: float = 1e-11             # d per pulse
    error_correction_efficiency: float = 1.1   # f, >= 1
    misalignment_error: float = 0.0            # e_a, [0, 0.5]
    visibility: float = 0.99                   # V, [0, 1] — NOT used by intensity_yield/qber_x_basis
```

Validated in `__post_init__`; raises `ValueError` on out-of-range inputs.

### `AnalyticalChannel`

```python
from snstfqkd.channel import AnalyticalChannel, ChannelParameters

channel = AnalyticalChannel(ChannelParameters(misalignment_error=0.2))
channel.intensity_yield(0.05, distance=400.0)   # S_mu, closed-form
channel.qber_x_basis(0.05, distance=400.0)      # E_mu^X, closed-form
```

Deterministic, instantaneous, no RNG. Use this whenever you need fast evaluation (e.g., inside an optimizer inner loop).

### `MonteCarloChannel`

```python
from snstfqkd.channel import MonteCarloChannel, ChannelParameters

channel = MonteCarloChannel(
    ChannelParameters(misalignment_error=0.2),
    mc_n_samples=2_000_000,     # trials per (mu, distance) evaluation
    mc_lambda_slice=1e-3,       # phase-slice half-width
    mc_seed=20260701,           # base seed for deterministic derivation
)
channel.intensity_yield(0.05, distance=400.0)   # S_mu, Monte-Carlo estimate
```

**Same `BaseChannel` interface as `AnalyticalChannel`** — anywhere you accept a `BaseChannel`, either works. `intensity_yield()` and `qber_x_basis()` are backed by a cached (`@lru_cache`) simulation; `z_window_yield()`, `gain()`, `vacuum_yield()`, `total_transmittance()`, `qber_z_basis()` use the *same closed-form formulas* as `AnalyticalChannel` (they were never patched to be MC-based upstream — don't "fix" this without checking `ARCHITECTURE.md` §13 first, it's intentional).

**Determinism:** calling `intensity_yield(mu, distance)` twice with the same arguments returns the identical result — the seed is derived from `(mc_seed, mu, distance, misalignment_error)`. This matters for `SNSOptimizer`'s `differential_evolution`, which would otherwise chase simulation noise.

Reduce `mc_n_samples` for fast iteration (e.g. `200_000` — used by the CLI's fast mode); use the default `2_000_000`+ for anything you'd publish.

### `channel.detector` — low-level MC engine

Only touch this if you're changing the Monte-Carlo physical model itself (rare — see §11.4). Key entry points:

```python
from snstfqkd.channel.detector import (
    mc_beamsplitter_xwindow_stats,   # production, cached, returns (S_mu, E_mu)
    mc_beamsplitter_diagnostics,     # diagnostic-only, uncached, returns a dict of intermediate quantities
    phase_slice_acceptance_fraction, # geometric acceptance rate for a given lambda_slice
    _derive_mc_seed,                 # deterministic seed derivation
)
```

`mc_beamsplitter_xwindow_stats` is `@lru_cache`d on its full positional-argument tuple — all arguments must be hashable (plain floats/ints, not numpy scalars) or the cache silently misses. Cast with `float(...)`/`int(...)` before calling if you're not sure.

---

## 4. Core API: `protocol/`

`SNSProtocolSimulator` is a **different abstraction** from `MonteCarloChannel` — it simulates explicit pulse-by-pulse trials rather than exposing the `BaseChannel` interface. Don't try to substitute one for the other.

```python
from snstfqkd.channel import AnalyticalChannel, ChannelParameters
from snstfqkd.protocol import SNSProtocolSimulator, evaluate_key_rate_from_simulation
from snstfqkd.security import DecoyStateAnalyzer, SNSKeyRateCalculator

channel = AnalyticalChannel(ChannelParameters(misalignment_error=0.15))
sim = SNSProtocolSimulator(
    channel,               # only .params and .total_transmittance() are used
    distance=300.0,
    epsilon=0.05,
    mu_signal=0.05,
    mu1=0.01,
    mu2=0.15,
    lambda_slice=1e-4,
    n_pulses_x=5_000_000,  # per-intensity trial count for X-windows
    n_pulses_z=5_000_000,  # trial count for Z-windows
    seed=42,
)

stats = sim.run()                      # -> ProtocolStatistics
# stats.S_mu0, stats.S_mu1, stats.S_mu2, stats.E_X_mu1, stats.S_Z

analyzer = DecoyStateAnalyzer()
calc = SNSKeyRateCalculator()
decoy, rate = evaluate_key_rate_from_simulation(stats, analyzer, calc)
```

You can also call the two halves separately for finer control:

```python
x1 = sim.simulate_x_windows(intensity=0.01)   # -> XWindowStats
z  = sim.simulate_z_windows()                  # -> ZWindowStats
```

`XWindowStats` exposes `S_mu`, `E_mu`, plus raw counts (`n_trials`, `n_phase_slice`, `n_effective`, `n_errors`, `n_double_click`, `n_no_click`) for diagnostics.

**Memory:** both `simulate_x_windows` and `simulate_z_windows` process pulses in fixed 2,000,000-pulse chunks internally, so `n_pulses_x=50_000_000` won't blow up memory — but it will take proportionally longer.

**Known discrepancy (not a bug):** `sim.simulate_z_windows().S_Z` differs from `channel.z_window_yield(mu_signal, distance)` — the protocol simulator captures the *true* ε-weighted mixture over all four Alice/Bob send/no-send combinations, while the analytical closed form implicitly assumes exactly one party always sends. See `protocol/protocol_simulator.py`'s module docstring.

### Pure helper functions

If you're writing a new protocol-level diagnostic, prefer composing these rather than re-deriving:

```python
from snstfqkd.protocol.phase_slice import sample_sliced_phases
from snstfqkd.protocol.events import bs_click_probs, sample_clicks, classify_x_window_events
from snstfqkd.protocol.pulse import sample_z_window_intensities, z_window_photon_number
```

All are plain numpy-in/numpy-out functions taking an explicit `np.random.Generator` — no hidden state, easy to unit test in isolation (see `tests/test_protocol_simulator.py::test_phase_slice_sampler_matches_rejection`).

---

## 5. Core API: `security/`

```python
from snstfqkd.security import DecoyStateAnalyzer, SNSKeyRateCalculator, binary_entropy, get_S_Z
```

### `DecoyStateAnalyzer` — Eq. (44), (45)

```python
analyzer = DecoyStateAnalyzer()
result = analyzer.three_intensity_decoy(mu1, mu2, s0, S_mu1, S_mu2, E_mu1)
# -> DecoyResult(s1, e1_ph, s0, valid, msg)

analyzer.validate_bounds(result)   # -> bool, sanity-checks s1/e1_ph ranges
```

Always check `result.valid` before using `result.s1`/`result.e1_ph` — invalid decoy configurations (e.g. `mu1 >= mu2`) return `valid=False` with `msg` explaining why, and `s1=0.0` as a safe default.

### `SNSKeyRateCalculator` — Eq. (4)

```python
calc = SNSKeyRateCalculator(f=1.1)   # f = error-correction efficiency

rate = calc.key_rate(epsilon, mu_signal, s1, e1_ph, S_Z, E_Z=0.0)

# Or sweep a distance array against a channel directly:
results = calc.key_rate_vs_distance(distances, channel, epsilon, mu_signal, mu1, mu2)
# -> List[KeyRateResult], each with .key_rate, .s1, .e1_ph, .S_Z, .distance
```

`key_rate_vs_distance` internally calls `get_S_Z(channel, mu_signal, d)` (see below) rather than `channel.intensity_yield()` — this is required, not optional; `intensity_yield` is the *X-window* observable and overestimates S_Z by roughly 2× if substituted.

### `get_S_Z` (from `security.bounds`)

```python
S_Z = get_S_Z(channel, mu_signal, distance)
```

Prefers `channel.z_window_yield()`; falls back to `channel.intensity_yield()` with a `UserWarning` if the channel doesn't implement `z_window_yield` (e.g. a minimal test double). If you're writing a new `BaseChannel` subclass, implement `z_window_yield()` — don't rely on the fallback.

---

## 6. Core API: `optimization/`

```python
from snstfqkd.optimization import SNSOptimizer, OptimizationResult

opt = SNSOptimizer(
    channel,
    f=1.1,
    epsilon_bounds=(1e-3, 0.5),
    mu_signal_bounds=(1e-3, 0.5),
    mu1_bounds=(1e-4, 0.3),
    mu2_bounds=(0.05, 0.5),
    de_maxiter=300,
)

result = opt.optimize_all(distance=300.0, method="differential_evolution", de_maxiter=200, seed=42)
# -> OptimizationResult(key_rate, epsilon, mu_signal, mu1, mu2, s1, e1_ph, method)
```

`method="grid"` is also supported (coarser, deterministic, no scipy dependency on the search itself — slower for the same resolution). Use `differential_evolution` unless you need reproducible-by-inspection grid search.

Single-parameter helpers, useful for diagnostics/plots rather than full joint optimization:

```python
opt.optimize_epsilon(distance, mu_signal, mu1, mu2, steps=50)
opt.optimize_mu_signal(distance, epsilon, mu1, mu2, steps=50)
opt.optimize_decoy_parameters(distance, epsilon, mu_signal, steps=30)  # -> (mu1, mu2)
```

**Cost model:** `optimize_all` with `differential_evolution` calls `_evaluate()` (one full `three_intensity_decoy` + `key_rate` computation) roughly `maxiter * popsize * 4` times. Against `AnalyticalChannel` this is cheap. Against `MonteCarloChannel`, each `_evaluate()` call triggers an MC simulation — cost scales with `mc_n_samples`. Keep `mc_n_samples` low (e.g. `200_000`) and/or `de_maxiter`/`de_popsize` low when optimizing against an MC channel, or expect long runtimes.

---

## 7. `config/` — Parameters & Presets

```python
from snstfqkd.config import make_channel_parameters, get_preset

# Build parameters with project defaults + overrides:
params = make_channel_parameters(misalignment_error=0.25)

# Or pull a named preset (returns (ChannelParameters, cfg_dict)):
params, cfg = get_preset("paper")
# cfg = {"channel_kwargs": {...}, "mu1": 0.01, "mu2": 0.15,
#        "mu_signal": 0.05, "epsilon": 0.05, "distances": [...]}
```

Available presets: `"paper"` (Wang 2018 defaults), `"low_misalignment"` (e_a=0), `"high_misalignment"` (e_a=0.35). Add new presets in `config/presets.py`'s `PRESETS` dict — see §11.1.

---

## 8. `visualization/` — Figure Functions

Every figure function follows the same shape: takes a channel (or channel class + params class for sweeps that need multiple channel instances), a distance/parameter array, optional keyword tunables, and `savepath: Optional[str]`. Returns `(fig, ax)` (or `(fig, ax, data...)` for comparison figures that also hand back the raw arrays for CSV export).

```python
import numpy as np
from snstfqkd.channel import AnalyticalChannel, ChannelParameters
from snstfqkd.visualization.analytical import figure_a8_s_mu1_vs_distance

channel = AnalyticalChannel(ChannelParameters(misalignment_error=0.15))
distances = np.linspace(0, 800, 100)

fig, ax = figure_a8_s_mu1_vs_distance(channel, distances, mu1=0.01, savepath="s_mu1.png")
```

If `savepath` is omitted, the figure is still returned (e.g. for interactive/notebook use) but nothing is written to disk.

Figures that need to instantiate *multiple* channels internally (e.g. A1/A2 sweep several misalignment values) take a **class + params class** instead of an instance:

```python
from snstfqkd.visualization.analytical import figure_a1_key_rate_vs_distance

figure_a1_key_rate_vs_distance(
    AnalyticalChannel, ChannelParameters, distances,
    misalignment_errors=(0.0, 0.15, 0.3), mu1=0.01, mu2=0.15,
    savepath="fig_a1.png",
)
```

Comparison figures (`visualization/comparison/`) take two already-constructed channel instances and return the raw data alongside the figure — this is how `workflows/comparison_workflow.py` builds its CSV:

```python
from snstfqkd.visualization.comparison import figure_c1_s_mu_comparison

fig, ax, s_mu_analytical, s_mu_mc = figure_c1_s_mu_comparison(
    analytical_channel, mc_channel, distances, mu1=0.01, savepath="c1.png"
)
```

Shared styling (`plt.rcParams`, color palette) lives in `visualization/_style.py`. Import it once at module top if you're adding a new figure module and want the same look:

```python
from snstfqkd.visualization._style import MISALIGNMENT_COLORS, LINESTYLES
```

---

## 9. `reports/` — CSV & Summary Reports

```python
from snstfqkd.reports import generate_distance_breakdown, write_breakdown_csv, write_summary_report

df = generate_distance_breakdown(
    channel, distances=[100, 300, 500, 700],
    optimize=True,             # per-distance optimize_all(), or...
    # optimize=False, fixed_params={"epsilon": 0.05, "mu_signal": 0.05, "mu1": 0.01, "mu2": 0.15},
    de_maxiter=100, de_popsize=8,
)
# df columns: distance_km, epsilon, mu_signal, mu1, mu2, gain, transmittance,
#             vacuum_yield, s0, S_mu1, S_mu2, E_mu1, s1, e1ph, S_Z,
#             positive_term, correction_term, R

write_breakdown_csv(df, "distance_breakdown.csv")
write_summary_report(df, "summary_report.txt", mode="analytical",
                      channel_params=channel.params, runtime_seconds=12.3)
```

`generate_distance_breakdown` works against **any** `BaseChannel` — this is what makes it usable for both the analytical and MC workflows without duplicating logic. `positive_term`/`correction_term`/`R` are the expanded terms of Eq. (4) (same formula as `SNSKeyRateCalculator.key_rate`, just with intermediates exposed for the CSV).

---

## 10. `workflows/` and `cli.py`

Each workflow function is a self-contained "generate everything for this mode" entry point, callable directly from Python (not just the CLI):

```python
from snstfqkd.workflows import run_analytical_workflow

df = run_analytical_workflow(
    output_root="my_results/analytical",
    preset="paper",
    fast=True,     # reduced optimizer/MC settings for speed
)
```

`fast=True` (the CLI default) trades accuracy for speed — fewer `differential_evolution` iterations, fewer MC samples, coarser sweeps. Use `fast=False` for anything you intend to publish or trust numerically.

The CLI (`cli.py`) is a thin `argparse` wrapper around these four functions — read it if you need to see exactly which workflow is invoked for which `--mode`, or to add a new mode (see §11.3).

---

## 11. Extension Recipes

### 11.1 Add a new named preset

Edit `config/presets.py`:

```python
PRESETS["long_haul"] = dict(
    channel_kwargs=dict(fiber_loss=0.18, misalignment_error=0.05),
    mu1=0.005, mu2=0.1, mu_signal=0.03, epsilon=0.03,
    distances=[200, 400, 600, 800, 1000, 1200, 1400],
)
```

No other changes needed — `get_preset("long_haul")` and `--preset long_haul` both pick it up immediately.

### 11.2 Add a new figure

1. Pick the right subpackage (`analytical/`, `mc/`, `protocol/`, or `comparison/`).
2. Write a function following the existing signature convention (channel(s) + array + kwargs + `savepath`), returning `(fig, ax)`.
3. Export it from that subpackage's `__init__.py`.
4. Wire it into the relevant `workflows/*_workflow.py` if it should run as part of `--mode <x>`.

```python
# visualization/analytical/my_new_figs.py
import matplotlib.pyplot as plt

def figure_a17_my_new_metric(channel, distances, savepath=None):
    vals = [channel.some_existing_method(d) for d in distances]
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(distances, vals)
    ax.set_xlabel("Distance (km)")
    ax.set_title("My New Metric")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax
```

**Only call existing `BaseChannel` methods or existing pure functions from `protocol.*`/`channel.detector`.** If the quantity you need genuinely doesn't exist anywhere yet, see §11.4 before adding new physics.

### 11.3 Add a new CLI mode / workflow

1. Write `workflows/my_workflow.py::run_my_workflow(output_root, preset, distances=None, fast=True) -> pd.DataFrame`, following the pattern in `workflows/mc_workflow.py` (create `figures/csv/reports` subdirs, generate figures, call `generate_distance_breakdown` + `write_breakdown_csv` + `write_summary_report`).
2. Export it from `workflows/__init__.py`.
3. In `cli.py`, add `"my_mode"` to the `--mode` choices and a corresponding `if args.mode in ("my_mode", "all"): run_my_workflow(...)` block.

### 11.4 Add a genuinely new physical quantity

This is the one case that touches the "no physics changed" constraint, so do it carefully:

- **If the quantity is derivable from existing observables** (e.g. a new figure of an existing `BaseChannel` method), just call the existing method — no new physics, see §11.2.
- **If it requires exposing an intermediate value an existing computation already computes internally but doesn't return** (this is what `mc_beamsplitter_diagnostics()` and `phase_slice_acceptance_fraction()` do), write a **new function that reproduces the existing math verbatim** and returns the additional field(s). Do not modify the production/cached function (`mc_beamsplitter_xwindow_stats`) — it's relied on by the optimizer's hot path and by the golden-value regression tests.
- **If it requires an actual new formula**, that's a physics change and needs sign-off + a new golden-value test, plus a paper/equation citation in the docstring. This should be rare — check `ARCHITECTURE.md` §2 (audit) and §5 (verification methodology) first to understand what's already been decided and why.

### 11.5 Add a new `BaseChannel` implementation

```python
from snstfqkd.channel.base import BaseChannel, ChannelParameters

class MyChannel(BaseChannel):
    def __init__(self, params: ChannelParameters):
        self.params = params

    def channel_transmittance(self, distance): ...
    def total_transmittance(self, distance): ...
    def vacuum_yield(self): ...
    def intensity_yield(self, intensity, distance): ...
    def z_window_yield(self, intensity, distance): ...
    def gain(self, intensity, distance): ...
    def qber_x_basis(self, intensity, distance): ...
    def qber_z_basis(self): ...
```

Implement all eight abstract methods (Python will raise `TypeError` at instantiation if you miss one). Once done, it's a drop-in replacement everywhere a `BaseChannel` is accepted: `security/`, `optimization/`, `visualization/analytical`+`mc`, `reports.generate_distance_breakdown`. Add golden-value tests analogous to `tests/test_analytical_channel.py` if this is meant to be a validated production channel.

---

## 12. Testing Conventions

```bash
pytest tests/ -q                          # everything (~2.5 min, dominated by MC/protocol tests)
pytest tests/test_analytical_channel.py   # fast (<1s), pure closed-form
pytest tests/test_mc_channel.py           # slower, does real MC sampling
pytest tests/test_cli.py                  # slowest, runs full "fast"-mode workflows end-to-end
```

**Golden-value tests** (`test_analytical_channel.py`, `test_mc_channel.py`, `test_optimizer.py`) hardcode numbers captured from the pre-refactor code. If you intentionally change a formula, you must:
1. Update the `GOLDEN` dict in the relevant test file.
2. Document the change and its physical justification in `ARCHITECTURE.md`.
3. Confirm the change doesn't silently affect downstream layers you didn't intend to touch (e.g. changing `AnalyticalChannel.intensity_yield` affects `security/`, `optimization/`, and most of `visualization/`).

If a golden-value test fails and you *didn't* intend to change physics, that's a regression — stop and diff against the original formula rather than updating the golden value.

**New pure functions** (anything in `protocol/phase_slice.py`, `protocol/events.py`, `protocol/pulse.py`, `channel/detector.py`) should get direct unit tests — they're plain numpy-in/numpy-out and cheap to test in isolation without spinning up a full channel or simulator.

**CLI/workflow tests** check file existence, not values — they exist to catch wiring breakage (missing import, wrong output path, etc.), not physics regressions. Don't add numeric assertions there; put those in the channel/security/optimizer test files instead.

---

## 13. Performance Notes & Gotchas

- **`MonteCarloChannel` inside an optimizer loop is expensive.** Each `differential_evolution` iteration calls `_evaluate()`, which calls `intensity_yield`/`qber_x_basis`, which runs a full MC simulation unless the `(mu, distance, e_a)` triple has been seen before (cache hit). Lower `mc_n_samples` for iteration; raise it only for the final run.
- **`lru_cache` on `mc_beamsplitter_xwindow_stats` requires hashable args.** Passing a numpy float64 instead of a Python float will still "work" (numpy floats are hashable) but can produce cache misses across what look like identical calls if the numpy dtype differs subtly. The channel classes already cast with `float(...)` before calling — keep doing this in any new caller.
- **`SNSProtocolSimulator` processes pulses in 2,000,000-pulse chunks** to bound peak memory. Don't remove this chunking when refactoring — `n_pulses_x=50_000_000` in one unchunked numpy array is ~400MB+ per intermediate array and will multiply across the half-dozen arrays in the hot loop.
- **`cli.py` calls `plt.close("all")` between workflow runs.** If you call workflow functions directly from a long-running Python process (not via the CLI) and generate many figures, do the same or you'll leak matplotlib figure objects (`RuntimeWarning: More than 20 figures have been opened`).
- **Fast vs. full CLI mode**: `--full` runs `MonteCarloChannel` at `2_000_000` samples and `SNSProtocolSimulator` at `5_000_000`+ pulses per point, across dozens of sweep points — expect this to take much longer than the default fast mode. Don't set `--full` as a default in CI.
- **`AnalyticalChannel` and `MonteCarloChannel` z_window_yield agree exactly** (same formula) — if you're debugging a mismatch in `z_window_yield`, the bug is not in the MC engine.
- **`SNSProtocolSimulator`'s S_Z will *not* match `channel.z_window_yield()`** — this is expected (see §4), don't "fix" it by making them agree without re-reading the module docstring first.

---

## 14. Coding Conventions for Contributors

- **Physics-touching code gets a docstring citation.** Any function implementing or approximating an equation from the paper should reference it (`Eq. (4)`, `Eq. (44)`, `Eq. (45)`, `Eq. (1)` for the phase-slice criterion, or a Wang et al. section/step number) directly in the docstring.
- **Prefer `BaseChannel` over concrete classes in type hints/signatures** unless you specifically need a Monte-Carlo-only or analytical-only feature.
- **New numeric constants go in `config/defaults.py`**, not as inline literals scattered across workflows/figures — this is precisely the anti-pattern the original codebase had (`breakdown.py`, `visualization.py`, `fig_01.py` each hardcoding their own defaults).
- **Pure numerical functions (no I/O, no plotting) belong in `channel/`, `protocol/`, or `security/`, not in `visualization/` or `workflows/`.** Visualization/workflow code should only call existing pure functions and format the output — this keeps the "expose a quantity without changing physics" pattern (§11.4) enforceable.
- **Every new `BaseChannel` method or `protocol.*` pure function should be independently unit-testable** without constructing a full workflow — if you find yourself needing a whole `SNSProtocolSimulator` instance to test a two-line function, extract the function (see how `phase_slice.py`/`events.py`/`pulse.py` were split out of `protocol_simulator.py` for the pattern to follow).
- **Don't add new module-level mutable state or singletons.** The `lru_cache` on `mc_beamsplitter_xwindow_stats` is deliberate (and its cache-key semantics are load-bearing for the optimizer); avoid introducing more caches without the same level of care around determinism.