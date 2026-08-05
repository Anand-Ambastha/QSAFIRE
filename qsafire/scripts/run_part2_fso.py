#!/usr/bin/env python3
"""Reproduce Part II of the notebook: the FSO + Gamma-Gamma atmospheric
turbulence extension of the SNS-TF-QKD engine (Parts 1-11).

Usage:
    python scripts/run_part2_fso.py [--no-plots]
"""
import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.experiments.fso_experiments import run_full_part2_pipeline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-plots", action="store_true", help="Skip matplotlib figure rendering.")
    args = parser.parse_args()

    run_full_part2_pipeline(show_plots=not args.no_plots)


if __name__ == "__main__":
    main()
