# snstfqkd

Research-grade SNS-TF-QKD simulation package: three independent, physics-unchanged
observable-generation layers (analytical closed-form, Monte-Carlo channel, and a
protocol-level pulse-by-pulse simulator), a shared decoy-state/key-rate security
layer (Eq. 44/45/4, unmodified), an optimizer, 46 figures across four categories,
CSV/report generation, and a unified CLI.

## Install

    pip install -e .

or with uv:

    uv pip install -e .

## Run

    uv run -m snstfqkd.main --mode analytical   # Figures A1-A16 + CSV + report
    uv run -m snstfqkd.main --mode mc            # Figures MC1-MC12 + CSV + report
    uv run -m snstfqkd.main --mode protocol      # Figures P1-P10 + CSV + report
    uv run -m snstfqkd.main --mode compare       # Figures C1-C8 + CSV + report
    uv run -m snstfqkd.main --mode all           # everything

By default the CLI runs in **fast** mode (reduced optimizer iterations / MC sample
counts) so a full `--mode all` run finishes in a couple of minutes. Pass `--full`
for publication-quality settings (slow).

Outputs are written to `results/<mode>/{figures,csv,reports}/`.

## Tests

    pytest tests/ -q

`tests/test_analytical_channel.py`, `tests/test_mc_channel.py`, and
`tests/test_optimizer.py` pin golden numerical values captured from the
**original, pre-refactor** code, so any accidental physics change will fail CI.

## Architecture

See the two audit documents delivered alongside this package for the full
rationale, dependency graph, and migration plan. In short:

    channel/       BaseChannel, AnalyticalChannel, MonteCarloChannel (Layers 1-2)
    protocol/      SNSProtocolSimulator, kept independent (Layer 3)
    security/      decoy_analysis.py (Eq. 44/45), key_rate.py (Eq. 4), bounds.py
    optimization/  SNSOptimizer
    visualization/ analytical/ mc/ protocol/ comparison/  (46 figures total)
    workflows/     one orchestration module per CLI mode
    reports/       distance_breakdown.csv + summary_report.txt generator
    config/        defaults.py, parameters.py, presets.py
    cli.py, main.py
