"""Top-level project configuration.

Physics constants live next to the equations that use them (``src/engine.py``,
``src/fso/config.py``) so that module reorganisation never risks silently
changing a numeric value used in a computation. This file only centralises
output-path configuration, shared across scripts.
"""

import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUTS_DIR = os.path.join(PROJECT_ROOT, "outputs")
FIGURES_DIR = os.path.join(OUTPUTS_DIR, "figures")
CSV_DIR = os.path.join(OUTPUTS_DIR, "csv")

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(CSV_DIR, exist_ok=True)
