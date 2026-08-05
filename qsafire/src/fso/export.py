"""Part 11 — Results Export.

CSV results and all figures (PNG + PDF) are written to ``outputs/csv/`` and
``outputs/figures/`` throughout the pipeline (see ``src/fso/metrics.py`` and
``src/fso/visualization.py``). This module writes the final short
human-readable summary report tying everything together, and lists exported
files.
"""

import glob
import pandas as pd

from src.engine import ETA_D, P_DARK, F_EC, DECOY_RATIO
from src.fso.config import (
    CSVDIR, FIGDIR, E_A, ALPHA_FSO_DB_KM, WAVELENGTH_M,
    TURBULENCE_CASES, DISTANCES_MAIN, DISTANCES_DET_EXTENDED,
)


def write_summary_report(summary_table, df_quad_val):
    report_lines = []
    report_lines.append("SNS-TF-QKD over FSO with Gamma-Gamma Atmospheric Turbulence")
    report_lines.append("=" * 60)
    report_lines.append(f"Generated: {pd.Timestamp.now().isoformat()}")
    report_lines.append("")
    report_lines.append("Fixed parameters:")
    report_lines.append(f"  Detector efficiency eta_d = {ETA_D}")
    report_lines.append(f"  Dark count probability p_dark = {P_DARK}")
    report_lines.append(f"  Error-correction factor f_EC = {F_EC}")
    report_lines.append(f"  Decoy ratio mu1/mu2 = {DECOY_RATIO}")
    report_lines.append(f"  Misalignment error e_a = {E_A}")
    report_lines.append(f"  Deterministic FSO attenuation alpha = {ALPHA_FSO_DB_KM} dB/km")
    report_lines.append(f"  Wavelength = {WAVELENGTH_M*1e9:.0f} nm")
    report_lines.append("")
    report_lines.append("Turbulence case studies (Cn2, m^-2/3):")
    for case, cn2 in TURBULENCE_CASES.items():
        report_lines.append(f"  {case:10s}: Cn2 = {cn2:.1e}")
    report_lines.append("")
    report_lines.append("Distance grids:")
    report_lines.append(f"  Main (realistic FSO)  : {DISTANCES_MAIN[0]:.2f} - {DISTANCES_MAIN[-1]:.2f} km one-way, {len(DISTANCES_MAIN)} points")
    report_lines.append(f"  Extended (det. only)  : {DISTANCES_DET_EXTENDED[0]:.2f} - {DISTANCES_DET_EXTENDED[-1]:.2f} km one-way, {len(DISTANCES_DET_EXTENDED)} points")
    report_lines.append("")
    report_lines.append("Maximum performance summary:")
    report_lines.append(summary_table.to_string(index=False))
    report_lines.append("")
    report_lines.append("Quadrature validation (Part 4): max relative error vs Monte Carlo across all")
    report_lines.append(f"tested distances and 20/40/80-point quadratures: {df_quad_val.rel_err.max():.4%} (target < 1%)")
    report_lines.append("")
    report_lines.append("Key qualitative findings:")
    report_lines.append(" 1. At the jointly-optimised operating point, R(eta) is locally CONVEX around")
    report_lines.append("    eta_det at short distances, so ensemble-averaged turbulence can transiently")
    report_lines.append("    RAISE the key rate above the deterministic reference (Jensen's inequality,")
    report_lines.append("    verified analytically and numerically - see Part 9 SKR figure discussion).")
    report_lines.append(" 2. The SNS protocol's key rate is extremely sensitive to the sending")
    report_lines.append("    probability epsilon; a naive choice (e.g. epsilon=0.05) yields NO secure")
    report_lines.append("    key at all, while the optimum sits in a narrow band around 0.02-0.03.")
    report_lines.append(" 3. Gamma-Gamma parameter alpha_g -> infinity, beta_g -> 1 as sigma_R^2 -> ")
    report_lines.append("    infinity (documented asymptotic feature of the Andrews-Phillips model);")
    report_lines.append("    handled numerically via a Gauss-Hermite fallback above shape parameter 150.")
    report_lines.append(" 4. Ensemble-optimal (eps*, mu'*) converge to the deterministic-optimal values")
    report_lines.append("    as distance grows and turbulence's relative influence on the optimum fades.")

    report_text = "\n".join(report_lines)
    with open(f"{CSVDIR}/summary_report.txt", "w") as f:
        f.write(report_text)
    print(report_text)
    return report_text


def list_exported_files():
    print("Exported CSV files:")
    for f in sorted(glob.glob(f"{CSVDIR}/*")):
        print(" ", f)
    print(f"\nExported figures ({len(glob.glob(f'{FIGDIR}/*.png'))} PNG + {len(glob.glob(f'{FIGDIR}/*.pdf'))} PDF):")
    for f in sorted(glob.glob(f"{FIGDIR}/*.png")):
        print(" ", f)
