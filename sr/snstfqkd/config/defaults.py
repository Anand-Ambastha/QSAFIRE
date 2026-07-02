# -*- coding: utf-8 -*-
"""
snstfqkd.config.defaults
===========================

Default numerical values used across workflows/CLI when the user does
not override them. These mirror the defaults already present in the
original ``ChannelParameters`` dataclass and in ``breakdown.py`` /
``fig_01.py`` / ``visualization.py`` — collected here so they live in
one place instead of being re-typed at each call site. No physics
changed.
"""

# Default channel parameters (Wang et al. 2018 paper values)
DEFAULT_FIBER_LOSS = 0.2                  # dB/km
DEFAULT_DETECTOR_EFFICIENCY = 0.8         # eta_d
DEFAULT_DARK_COUNT_RATE = 1e-11           # d per pulse
DEFAULT_ERROR_CORRECTION_EFFICIENCY = 1.1  # f
DEFAULT_MISALIGNMENT_ERROR = 0.15         # e_a (breakdown.py default)
DEFAULT_VISIBILITY = 0.99                 # V

# Default decoy intensities
DEFAULT_MU1 = 0.01
DEFAULT_MU2 = 0.15
DEFAULT_MU_SIGNAL = 0.05
DEFAULT_EPSILON = 0.05

# Default distance sweep (breakdown.py / visualization.py)
DEFAULT_DISTANCES_KM = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1300]
DEFAULT_MAX_DISTANCE = 800.0
DEFAULT_DISTANCE_POINTS = 100

# Default misalignment sweep (visualization.py Fig. 2)
DEFAULT_MISALIGNMENTS = (0.00, 0.15, 0.25, 0.35, 0.45)

# Default optimizer settings
DEFAULT_DE_MAXITER = 200
DEFAULT_DE_POPSIZE = 12
DEFAULT_SEED = 42

# Default Monte-Carlo channel settings (sns_tfqkd_mc.py defaults)
DEFAULT_MC_N_SAMPLES = 2_000_000
DEFAULT_MC_LAMBDA_SLICE = 1.0e-3
DEFAULT_MC_SEED = 20260701

# Default protocol-simulator settings (protocol_simulator.py defaults)
DEFAULT_PROTOCOL_LAMBDA_SLICE = 1e-4
DEFAULT_PROTOCOL_N_PULSES_X = 50_000_000
DEFAULT_PROTOCOL_N_PULSES_Z = 50_000_000

# fig_01.py hybrid MC/analytical crossover distance
DEFAULT_MC_CROSSOVER_KM = 200.0
