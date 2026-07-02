"""
test_cli.py

End-to-end smoke test of the CLI: verifies each mode produces the
expected figures/CSV/report files. Uses --output-root pointed at a
pytest tmp_path and 'fast' (default) settings for speed. Verifies
wiring/file output only — does not re-derive physics values (see
test_analytical_channel.py / test_mc_channel.py / test_optimizer.py /
test_security.py for those).
"""
import os
import pytest

from snstfqkd.cli import main


def test_cli_analytical_mode(tmp_path):
    out = str(tmp_path / "results")
    main(["--mode", "analytical", "--output-root", out])
    assert os.path.exists(os.path.join(out, "analytical", "csv", "distance_breakdown.csv"))
    assert os.path.exists(os.path.join(out, "analytical", "reports", "summary_report.txt"))
    assert os.path.exists(os.path.join(out, "analytical", "figures", "fig_a1_key_rate_vs_distance.png"))


def test_cli_mc_mode(tmp_path):
    out = str(tmp_path / "results")
    main(["--mode", "mc", "--output-root", out])
    assert os.path.exists(os.path.join(out, "mc", "csv", "distance_breakdown.csv"))
    assert os.path.exists(os.path.join(out, "mc", "reports", "summary_report.txt"))
    assert os.path.exists(os.path.join(out, "mc", "figures", "fig_mc1_key_rate_vs_distance.png"))


def test_cli_protocol_mode(tmp_path):
    out = str(tmp_path / "results")
    main(["--mode", "protocol", "--output-root", out])
    assert os.path.exists(os.path.join(out, "protocol", "csv", "distance_breakdown.csv"))
    assert os.path.exists(os.path.join(out, "protocol", "reports", "summary_report.txt"))
    assert os.path.exists(os.path.join(out, "protocol", "figures", "fig_p1_z_window_fraction.png"))


def test_cli_compare_mode(tmp_path):
    out = str(tmp_path / "results")
    main(["--mode", "compare", "--output-root", out])
    assert os.path.exists(os.path.join(out, "comparison", "csv", "distance_breakdown.csv"))
    assert os.path.exists(os.path.join(out, "comparison", "reports", "summary_report.txt"))
    assert os.path.exists(os.path.join(out, "comparison", "figures", "fig_c1_s_mu_comparison.png"))
