#!/usr/bin/env python3
"""Top-level reproduction entry point.

Running this script reproduces every currently-migrated part of the source
notebook end-to-end, without Jupyter:

  Part I  - fibre SNS-TF-QKD simulation (Sections 1-12)
  Part II - FSO + Gamma-Gamma atmospheric turbulence extension (Parts 1-11)

NOTE: the source notebook also contains Parts III-XIII (satellite orbital
mechanics, atmospheric transmittance, slant-path turbulence, seasonal
extension, PAT/optics losses, weather-gated availability, joint
optimisation, multi-city links, and per-link asymmetric metrics). Those
parts are not yet migrated into this package; see docs/module_documentation.md
for migration status.

Usage:
    python scripts/train.py [--no-plots]
"""
import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import scripts.run_part1_fiber as run_part1_fiber  # noqa: E402
from src.experiments.fso_experiments import run_full_part2_pipeline  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-plots", action="store_true", help="Skip matplotlib figure rendering.")
    args = parser.parse_args()
    show_plots = not args.no_plots

    print("\n" + "#" * 70)
    print("# PART I - Fibre SNS-TF-QKD reproduction")
    print("#" * 70)
    run_part1_fiber.main()

    print("\n" + "#" * 70)
    print("# PART II - FSO + Gamma-Gamma turbulence extension")
    print("#" * 70)
    run_full_part2_pipeline(show_plots=show_plots)


if __name__ == "__main__":
    main()
