"""Shared configuration/constants for the Part II FSO + Gamma-Gamma turbulence
extension (notebook Section 13 imports and Part 1/2 constants).
"""

import os
import numpy as np

# ---------------------------------------------------------------------------
# Section 13: output directories
# ---------------------------------------------------------------------------
FIGDIR = "outputs/figures"
CSVDIR = "outputs/csv"
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(CSVDIR, exist_ok=True)
CASE_COLORS = {"weak": "#2a9d8f", "moderate": "#e76f51", "strong": "#8e44ad"}

# ---------------------------------------------------------------------------
# Part 1: Model A (deterministic FSO) constant
# ---------------------------------------------------------------------------
ALPHA_FSO_DB_KM = 0.43

# ---------------------------------------------------------------------------
# Part 1: Model B (Gamma-Gamma turbulence) constants
# ---------------------------------------------------------------------------
WAVELENGTH_M = 1550e-9
K_WAVENUMBER = 2 * np.pi / WAVELENGTH_M

CN2_WEAK = 5e-15       # m^-2/3, weak turbulence case study
CN2_MODERATE = 1e-14   # m^-2/3, moderate turbulence case study
CN2_STRONG = 5e-14     # m^-2/3, strong turbulence case study
TURBULENCE_CASES = {"weak": CN2_WEAK, "moderate": CN2_MODERATE, "strong": CN2_STRONG}

# ---------------------------------------------------------------------------
# Part 2: distance grids and fixed misalignment error
# ---------------------------------------------------------------------------
DISTANCES_MAIN = np.round(np.geomspace(0.3, 200, 16), 3)
DISTANCES_DET_EXTENDED = np.round(np.geomspace(0.3, 250, 22), 2)
E_A = 0.15

# ---------------------------------------------------------------------------
# Parts 5-7: joint optimisation bounds
# ---------------------------------------------------------------------------
EPS_BOUNDS = (0.01, 0.99)
MU_BOUNDS = (0.01, 1.0)

# ---------------------------------------------------------------------------
# Part 3: Monte Carlo default sample size
# ---------------------------------------------------------------------------
N_MC_DEFAULT = 100_000

# ---------------------------------------------------------------------------
# Part 10: secure-key-rate threshold
# ---------------------------------------------------------------------------
SECURE_THRESHOLD = 1e-12

# ---------------------------------------------------------------------------
# Part 6: naive epsilon baselines for the improvement-factor figure
# ---------------------------------------------------------------------------
EPS_NAIVE_FAIL, EPS_NAIVE_MILD = 0.05, 0.01
