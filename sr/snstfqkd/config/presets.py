# -*- coding: utf-8 -*-
"""
snstfqkd.config.presets
==========================

Named parameter presets used by the CLI/workflows. "paper" reproduces
the Wang et al. 2018 default operating point used throughout the
original codebase (breakdown.py, visualization.py, fig_01.py).
"""

from snstfqkd.config.parameters import make_channel_parameters
from snstfqkd.config import defaults as D

PRESETS = {
    "paper": dict(
        channel_kwargs=dict(),  # pure defaults
        mu1=D.DEFAULT_MU1,
        mu2=D.DEFAULT_MU2,
        mu_signal=D.DEFAULT_MU_SIGNAL,
        epsilon=D.DEFAULT_EPSILON,
        distances=D.DEFAULT_DISTANCES_KM,
    ),
    "low_misalignment": dict(
        channel_kwargs=dict(misalignment_error=0.0),
        mu1=D.DEFAULT_MU1,
        mu2=D.DEFAULT_MU2,
        mu_signal=D.DEFAULT_MU_SIGNAL,
        epsilon=D.DEFAULT_EPSILON,
        distances=D.DEFAULT_DISTANCES_KM,
    ),
    "high_misalignment": dict(
        channel_kwargs=dict(misalignment_error=0.35),
        mu1=D.DEFAULT_MU1,
        mu2=D.DEFAULT_MU2,
        mu_signal=D.DEFAULT_MU_SIGNAL,
        epsilon=D.DEFAULT_EPSILON,
        distances=D.DEFAULT_DISTANCES_KM,
    ),
}


def get_preset(name: str = "paper"):
    if name not in PRESETS:
        raise KeyError(f"Unknown preset '{name}'. Available: {list(PRESETS)}")
    cfg = PRESETS[name]
    params = make_channel_parameters(**cfg["channel_kwargs"])
    return params, cfg
