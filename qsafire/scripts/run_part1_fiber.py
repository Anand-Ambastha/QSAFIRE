#!/usr/bin/env python3
"""Reproduce Part I of the notebook: the fibre-based SNS-TF-QKD simulation
(Sections 1-12), matching Wang, Yu & Hu (2018), Section III.

Usage:
    python scripts/run_part1_fiber.py [--no-plots]
"""
import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.experiments.fiber_experiments import (
    demo_channel_scaling, demo_gain_and_error, demo_z_qber,
    validate_s1_decoy_estimate, demo_e1ph_recovery,
    quick_check_key_rate, quick_check_optimizer,
    figure1_sweep, figure2_sweep, validation_tables,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-plots", action="store_true", help="Skip matplotlib figure rendering.")
    args = parser.parse_args()
    show_plots = not args.no_plots

    print("=== Section 3: channel scaling ===")
    print(demo_channel_scaling())

    print("\n=== Section 4: gain/error illustration ===")
    demo_gain_and_error(show_plots)

    print("\n=== Section 5: Z-basis QBER illustration ===")
    demo_z_qber(show_plots)

    print("\n=== Section 6: decoy-state s1 validation ===")
    print(validate_s1_decoy_estimate())

    print("\n=== Section 7: e1^ph recovery illustration ===")
    print(demo_e1ph_recovery())

    print("\n=== Section 8: key-rate quick check ===")
    quick_check_key_rate()

    print("\n=== Section 9: optimiser quick check ===")
    quick_check_optimizer()

    print("\n=== Section 10: Figure 1 reproduction ===")
    figure1_sweep(show_plots)

    print("\n=== Section 11: Figure 2 reproduction ===")
    figure2_sweep(show_plot=show_plots)

    print("\n=== Section 12: validation tables ===")
    df_validation, df_focus = validation_tables()
    print(df_validation)


if __name__ == "__main__":
    main()
