# -*- coding: utf-8 -*-
"""Shared matplotlib style, moved verbatim from the original visualization.py."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.linewidth": 1.2,
    "lines.linewidth": 1.8,
    "legend.fontsize": 9,
    "legend.framealpha": 0.85,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

MISALIGNMENT_COLORS = ["#1a1a2e", "#16213e", "#0f3460", "#533483", "#e94560"]
LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]
