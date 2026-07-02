# -*- coding: utf-8 -*-
"""
snstfqkd.security.bounds
===========================

Shared helper used by both ``key_rate.py`` and ``optimizer.py`` in the
original codebase — ``_get_S_Z`` was duplicated identically in both
files. Deduplicated here with no change to the logic.
"""

import warnings


def get_S_Z(channel, mu_signal: float, distance: float) -> float:
    """
    Return the Z-window counting rate S_Z.

    Uses channel.z_window_yield() if available (patched ChannelModel),
    otherwise falls back to channel.intensity_yield() with a warning.

    Wang 2018 Eq. (4): S_Z is the Z-window single-click rate, which uses
    one-arm BS statistics (one party sends, the other sends vacuum).
    This is different from the X-window yield returned by intensity_yield().
    """
    if hasattr(channel, 'z_window_yield'):
        return channel.z_window_yield(mu_signal, distance)
    else:
        warnings.warn(
            "channel.z_window_yield() not found; falling back to "
            "intensity_yield() for S_Z. This overestimates S_Z.",
            stacklevel=3,
        )
        return channel.intensity_yield(mu_signal, distance)
