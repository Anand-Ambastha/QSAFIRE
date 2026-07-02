# -*- coding: utf-8 -*-
"""
snstfqkd.cli
==============

Unified command-line interface.

    uv run -m snstfqkd.main --mode analytical
    uv run -m snstfqkd.main --mode mc
    uv run -m snstfqkd.main --mode protocol
    uv run -m snstfqkd.main --mode compare
    uv run -m snstfqkd.main --mode all
"""

import argparse

import matplotlib.pyplot as plt

from snstfqkd.workflows import (
    run_analytical_workflow, run_mc_workflow,
    run_protocol_workflow, run_comparison_workflow,
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="snstfqkd",
        description="SNS-TF-QKD research-grade simulation & figure generation CLI.",
    )
    p.add_argument(
        "--mode", choices=["analytical", "mc", "protocol", "compare", "all"],
        default="all", help="Which workflow(s) to run.",
    )
    p.add_argument("--preset", default="paper", help="Named parameter preset (see config.presets).")
    p.add_argument("--output-root", default="results", help="Root output directory.")
    p.add_argument(
        "--full", action="store_true",
        help="Disable 'fast' mode: use full optimizer iterations / MC sample counts "
             "(slow, publication-quality). Default is fast mode for quick CLI runs.",
    )
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    fast = not args.full

    ran = []
    if args.mode in ("analytical", "all"):
        run_analytical_workflow(output_root=f"{args.output_root}/analytical",
                                 preset=args.preset, fast=fast)
        plt.close("all")
        ran.append("analytical")
    if args.mode in ("mc", "all"):
        run_mc_workflow(output_root=f"{args.output_root}/mc", preset=args.preset, fast=fast)
        plt.close("all")
        ran.append("mc")
    if args.mode in ("protocol", "all"):
        run_protocol_workflow(output_root=f"{args.output_root}/protocol",
                               preset=args.preset, fast=fast)
        plt.close("all")
        ran.append("protocol")
    if args.mode in ("compare", "all"):
        run_comparison_workflow(output_root=f"{args.output_root}/comparison",
                                 preset=args.preset, fast=fast)
        plt.close("all")
        ran.append("compare")

    print(f"snstfqkd: completed workflows: {', '.join(ran)}")
    print(f"Outputs written under: {args.output_root}/")


if __name__ == "__main__":
    main()
