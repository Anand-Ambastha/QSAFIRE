# -*- coding: utf-8 -*-
"""
Figure MC12 — Runtime vs Number of MC Samples

Pure performance diagnostic: times ``mc_beamsplitter_xwindow_stats``
(the production MC engine, unmodified) at a fixed operating point over
a range of ``n_samples`` values. No physics involved.
"""

import time
from typing import Optional, Sequence

import numpy as np
import matplotlib.pyplot as plt

from snstfqkd.channel.detector import mc_beamsplitter_xwindow_stats


def figure_mc12_runtime_vs_samples(
    mc_channel,
    sample_counts: Sequence[int] = (10_000, 50_000, 200_000, 500_000, 1_000_000, 2_000_000),
    distance: float = 300.0,
    intensity: float = 0.05,
    savepath: Optional[str] = None,
):
    runtimes = []
    for n in sample_counts:
        mc_beamsplitter_xwindow_stats.cache_clear()
        t0 = time.perf_counter()
        mc_beamsplitter_xwindow_stats(
            mc_channel.params.fiber_loss,
            mc_channel.params.detector_efficiency,
            mc_channel.params.dark_count_rate,
            mc_channel.params.misalignment_error,
            float(intensity), float(distance),
            int(n), mc_channel._mc_lambda_slice, mc_channel._mc_seed,
        )
        runtimes.append(time.perf_counter() - t0)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.loglog(sample_counts, runtimes, "o-", ms=4, color="tab:brown")
    ax.set_xlabel("Number of MC samples")
    ax.set_ylabel("Runtime (s)")
    ax.set_title("Runtime vs Number of MC Samples")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")
    plt.tight_layout()
    if savepath:
        fig.savefig(savepath)
    return fig, ax
