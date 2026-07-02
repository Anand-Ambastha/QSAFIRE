# -*- coding: utf-8 -*-
"""
snstfqkd.config.parameters
=============================

Convenience factory for building ``ChannelParameters`` /
``SimulationParameters`` from the defaults in ``config.defaults``,
with keyword overrides. No new physics — just a constructor
convenience layer used by workflows and the CLI.
"""

from snstfqkd.channel.base import ChannelParameters, SimulationParameters
from snstfqkd.config import defaults as D


def make_channel_parameters(**overrides) -> ChannelParameters:
    """Build ChannelParameters from project defaults, with overrides."""
    base = dict(
        fiber_loss=D.DEFAULT_FIBER_LOSS,
        detector_efficiency=D.DEFAULT_DETECTOR_EFFICIENCY,
        dark_count_rate=D.DEFAULT_DARK_COUNT_RATE,
        error_correction_efficiency=D.DEFAULT_ERROR_CORRECTION_EFFICIENCY,
        misalignment_error=D.DEFAULT_MISALIGNMENT_ERROR,
        visibility=D.DEFAULT_VISIBILITY,
    )
    base.update(overrides)
    return ChannelParameters(**base)


def make_simulation_parameters(**overrides) -> SimulationParameters:
    """Build SimulationParameters from project defaults, with overrides."""
    base = dict(
        min_distance=0.0,
        max_distance=D.DEFAULT_MAX_DISTANCE,
        distance_points=D.DEFAULT_DISTANCE_POINTS,
    )
    base.update(overrides)
    return SimulationParameters(**base)
