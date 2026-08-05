"""Part II orchestrator: reproduces the FSO/turbulence extension notebook
cells in order (Parts 1-11), producing the same tables, validation checks,
figures, and exported files.
"""

import numpy as np
import pandas as pd
from scipy.integrate import quad

from src.fso.config import (
    DISTANCES_MAIN, DISTANCES_DET_EXTENDED, E_A, TURBULENCE_CASES,
    CN2_MODERATE, N_MC_DEFAULT,
)
from src.fso.channel_models import (
    gamma_gamma_params, gamma_gamma_pdf, sample_gamma_gamma, rytov_variance,
)
from src.fso.monte_carlo import monte_carlo_convergence
from src.fso.quadrature import gauss_laguerre_ensemble_rate, quadrature_validation
from src.fso.optimization import instantaneous_maps
from src.fso.metrics import build_deterministic_table, build_gamma_gamma_table, maximum_performance_summary
from src.fso.visualization import (
    plot_rytov_variance_vs_distance, plot_turbulence_regime_map,
    plot_all_deterministic_figures, plot_all_gamma_gamma_figures,
    plot_optimal_param_vs_distance, plot_epsilon_improvement,
    plot_instantaneous_map, plot_joint_optimization_landscape,
    plot_quadrature_validation, plot_monte_carlo_convergence,
)
from src.fso.export import write_summary_report, list_exported_files


def validate_gamma_gamma_sampler():
    """Cells 33-34: validate the Gamma-Gamma sampler against theoretical
    variance and PDF normalisation."""
    rng_check = np.random.default_rng(42)
    rows = []
    for s2 in [0.1, 0.5, 1.0, 3.0, 8.0, 15.0]:
        ag, bg = gamma_gamma_params(s2)
        h = sample_gamma_gamma(ag, bg, 1_000_000, rng_check)
        var_theory = 1 / ag + 1 / bg + 1 / (ag * bg)
        rows.append((s2, ag, bg, h.mean(), h.var(), var_theory))
    df_gg_check = pd.DataFrame(rows, columns=["sigma_R2", "alpha_g", "beta_g", "mean(h) [MC]", "var(h) [MC]", "var(h) [theory]"])
    print(df_gg_check)

    ag, bg = gamma_gamma_params(3.0)
    pdf_integral, _ = quad(lambda h: gamma_gamma_pdf(h, ag, bg), 1e-8, 200, limit=200)
    pdf_mean, _ = quad(lambda h: h * gamma_gamma_pdf(h, ag, bg), 1e-8, 200, limit=200)
    print(f"PDF integral (want 1): {pdf_integral:.6f}   PDF mean (want 1): {pdf_mean:.6f}")
    print("\nNote: alpha_g -> infinity and beta_g -> 1 as sigma_R^2 -> infinity is the known")
    print("asymptotic (fully-saturated / negative-exponential) limit of the Gamma-Gamma model;")
    print("Var(h) is NOT monotonic in sigma_R^2 (it peaks around sigma_R^2~1, dips through the")
    print("moderate/strong 'focusing' regime, then re-approaches 1 only in extreme saturation")
    print("far beyond typical terrestrial FSO ranges) - this is a genuine, documented feature")
    print("of the standard Andrews-Phillips formulas, not an artifact of this implementation.")
    return df_gg_check


def part1_2_channel_and_regime(show_plots=True):
    """Cells 32-37: channel models + Gamma-Gamma sampler validation + regime
    map."""
    validate_gamma_gamma_sampler()
    if show_plots:
        plot_rytov_variance_vs_distance()
        plot_turbulence_regime_map()


def part3_4_mc_and_quadrature_validation(show_plots=True):
    """Cells 40-43: quadrature-vs-Monte-Carlo validation and MC convergence
    diagnostic."""
    val_rows = []
    for d in DISTANCES_MAIN:
        s2 = rytov_variance(d, CN2_MODERATE)
        mc, gl_rows = quadrature_validation(0.02, 0.3, d, E_A, s2, n_mc=150_000, n_points_list=(20, 40, 80))
        for r in gl_rows:
            val_rows.append(dict(d_km=d, sigma_R2=s2, **r))
    df_quad_val = pd.DataFrame(val_rows)
    print(df_quad_val.pivot(index="d_km", columns="n_points", values="rel_err"))

    if show_plots:
        plot_quadrature_validation(df_quad_val)
    print(f"Max relative error across all tested points/distances: {df_quad_val.rel_err.max():.4%}  (target < 1%)")

    # Monte Carlo convergence diagnostic (cell 43)
    n_list = np.unique(np.logspace(1, 5, 30).astype(int)).tolist()
    d_demo, s2_demo = 6.0, rytov_variance(6.0, CN2_MODERATE)
    n_arr, running_mean = monte_carlo_convergence(0.02, 0.3, d_demo, E_A, s2_demo, n_list, seed=0)
    true_val = gauss_laguerre_ensemble_rate(0.02, 0.3, d_demo, E_A, s2_demo, n_points=60)
    if show_plots:
        plot_monte_carlo_convergence(n_arr, running_mean, true_val, d_demo)

    return df_quad_val


def part5_11_full_pipeline(show_plots=True):
    """Cells 48-77: metrics tables, all Part 9 figures, Part 10 summary, and
    Part 11 export."""
    df_det = build_deterministic_table(DISTANCES_MAIN)
    df_det_ext = build_deterministic_table(DISTANCES_DET_EXTENDED)  # Part 10 secure-distance only
    df_gg = build_gamma_gamma_table(DISTANCES_MAIN, n_mc=100_000)

    df_det.to_csv("outputs/csv/deterministic_fso_results.csv", index=False)
    df_gg.to_csv("outputs/csv/gamma_gamma_results.csv", index=False)

    print("Deterministic FSO table:")
    print(df_det.round(6))
    print("Gamma-Gamma table (head):")
    print(df_gg.round(6).head(10))

    if show_plots:
        plot_all_deterministic_figures(df_det)
        plot_all_gamma_gamma_figures(df_det, df_gg)
        plot_optimal_param_vs_distance(df_det, df_gg, "mu", "Optimal $\\mu'$",
                                        "Optimal signal intensity vs distance", "alpha_optimization_curve")
        plot_optimal_param_vs_distance(df_det, df_gg, "eps", "Optimal $\\epsilon$",
                                        "Optimal sending probability vs distance", "epsilon_optimization_curve")
        plot_epsilon_improvement(df_gg)

        df_inst = pd.DataFrame(instantaneous_maps(DISTANCES_MAIN, E_A, lambda d: rytov_variance(d, CN2_MODERATE)))
        df_inst.to_csv("outputs/csv/instantaneous_optimization_map.csv", index=False)
        plot_instantaneous_map(df_inst, "mu_inst", "Optimal $\\mu'_{inst}$",
                                "Instantaneous alpha map (moderate turbulence)", "instantaneous_alpha_map")
        plot_instantaneous_map(df_inst, "eps_inst", "Optimal $\\epsilon_{inst}$",
                                "Instantaneous epsilon map (moderate turbulence)", "instantaneous_epsilon_map")

        plot_joint_optimization_landscape(d_demo=6.0)

    summary_table = maximum_performance_summary(df_det_ext, df_gg)
    summary_table.to_csv("outputs/csv/max_performance_summary.csv", index=False)
    print(summary_table)

    return df_det, df_det_ext, df_gg, summary_table


def run_full_part2_pipeline(show_plots=True):
    """Runs the entire Part II FSO extension end-to-end, in notebook order,
    and writes the final summary report + file listing."""
    part1_2_channel_and_regime(show_plots)
    df_quad_val = part3_4_mc_and_quadrature_validation(show_plots)
    df_det, df_det_ext, df_gg, summary_table = part5_11_full_pipeline(show_plots)
    write_summary_report(summary_table, df_quad_val)
    list_exported_files()
    return dict(df_det=df_det, df_det_ext=df_det_ext, df_gg=df_gg,
                summary_table=summary_table, df_quad_val=df_quad_val)
