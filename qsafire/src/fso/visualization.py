"""Part 2 (regime-map figures) and Part 9 — Publication-Quality Figures.

Every figure below preserves the exact colors, labels, legends, scales,
annotations, titles and layout of the corresponding notebook cell. Figures
are saved as PNG + PDF into ``FIGDIR`` exactly as in the notebook.
"""

import numpy as np
import matplotlib.pyplot as plt

from src.fso.config import FIGDIR, CASE_COLORS, TURBULENCE_CASES, DISTANCES_MAIN, E_A
from src.fso.channel_models import (
    rytov_variance, gamma_gamma_params, classify_turbulence_regime,
    channel_loss_dB, sample_gamma_gamma,
)
from src.fso.quadrature import gauss_laguerre_ensemble_rate
from src.fso.optimization import joint_optimize_ensemble, EPS_BOUNDS, MU_BOUNDS


def save_fig(fig, name):
    fig.savefig(f"{FIGDIR}/{name}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{FIGDIR}/{name}.pdf", bbox_inches='tight')


# ---------------------------------------------------------------------------
# Part 2 figures
# ---------------------------------------------------------------------------
def plot_rytov_variance_vs_distance(distances_main=DISTANCES_MAIN):
    fig, ax = plt.subplots(figsize=(7, 5))
    for case, color in CASE_COLORS.items():
        s2 = rytov_variance(distances_main, TURBULENCE_CASES[case])
        ax.plot(distances_main, s2, 'o-', color=color, ms=4, label=f"{case} ($C_n^2$={TURBULENCE_CASES[case]:.0e})")
    ax.axhline(1.0, color='gray', ls=':', lw=1)
    ax.axhline(10.0, color='gray', ls=':', lw=1)
    ax.text(distances_main[0], 1.15, "weak / moderate boundary", fontsize=8, color='gray')
    ax.text(distances_main[0], 11.5, "moderate / strong boundary", fontsize=8, color='gray')
    ax.set_yscale('log')
    ax.set_xlabel("One-way FSO distance $d$ (km)")
    ax.set_ylabel(r"Rytov variance $\sigma_R^2$")
    ax.set_title("Rytov variance vs distance")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "rytov_variance_vs_distance")
    plt.show()
    return fig


def plot_turbulence_regime_map(distances_main=DISTANCES_MAIN):
    regime_code = {"weak": 0, "moderate": 1, "strong": 2}
    cmap = plt.matplotlib.colors.ListedColormap(["#2a9d8f", "#e9c46a", "#e76f51"])
    cases = list(TURBULENCE_CASES.keys())
    grid = np.zeros((len(cases), len(distances_main)))
    for i, case in enumerate(cases):
        s2 = rytov_variance(distances_main, TURBULENCE_CASES[case])
        grid[i, :] = [regime_code[r] for r in classify_turbulence_regime(s2)]

    fig, ax = plt.subplots(figsize=(8, 3.2))
    im = ax.imshow(grid, aspect='auto', cmap=cmap, vmin=-0.5, vmax=2.5,
                    extent=[distances_main[0], distances_main[-1], -0.5, len(cases) - 0.5], origin='lower')
    ax.set_yticks(range(len(cases)))
    ax.set_yticklabels(cases)
    ax.set_xlabel("One-way FSO distance $d$ (km)")
    ax.set_title("Turbulence regime map (Gamma-Gamma validity classification)")
    cbar = fig.colorbar(im, ax=ax, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(["weak", "moderate", "strong/saturation"])
    fig.tight_layout()
    save_fig(fig, "turbulence_regime_map")
    plt.show()
    return fig


# ---------------------------------------------------------------------------
# Part 9: Deterministic FSO figures
# ---------------------------------------------------------------------------
def plot_det_vs_distance(df_det, ycol, ylabel, title, name, logy=True):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(df_det["d_km"], df_det[ycol], 'o-', color='#1d3557', ms=4)
    if logy:
        ax.set_yscale('log')
    ax.set_xlabel("One-way FSO distance $d$ (km)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, name)
    plt.show()
    return fig


def plot_all_deterministic_figures(df_det):
    plot_det_vs_distance(df_det, "R", "Key rate $R$ (bits/pulse)", "Deterministic FSO: SKR vs distance", "det_skr_vs_distance")
    plot_det_vs_distance(df_det, "Qxx", "$Q_{xx}$ (gain)", "Deterministic FSO: $Q_{xx}$ vs distance", "det_qxx_vs_distance", logy=False)
    plot_det_vs_distance(df_det, "Exx", "$E_{xx}$ (QBER)", "Deterministic FSO: $E_{xx}$ vs distance", "det_exx_vs_distance", logy=False)
    plot_det_vs_distance(df_det, "s1", "$s_1$ (single-photon yield)", "Deterministic FSO: $s_1$ vs distance", "det_s1_vs_distance", logy=False)
    plot_det_vs_distance(df_det, "loss_dB", "Channel loss (dB)", "Deterministic FSO: channel loss vs distance", "det_loss_vs_distance", logy=False)
    plot_det_vs_distance(df_det, "mu", "Optimal $\\mu'$", "Deterministic FSO: optimal signal intensity vs distance", "det_mu_vs_distance", logy=False)
    plot_det_vs_distance(df_det, "eps", "Optimal $\\epsilon$", "Deterministic FSO: optimal sending probability vs distance", "det_eps_vs_distance", logy=False)


# ---------------------------------------------------------------------------
# Part 9: Gamma-Gamma figures
# ---------------------------------------------------------------------------
def plot_gg_vs_distance(df_gg, ycol, ylabel, title, name, logy=True, det_reference=None, ci=False):
    fig, ax = plt.subplots(figsize=(7, 5))
    for case, color in CASE_COLORS.items():
        sub = df_gg[df_gg["turbulence_case"] == case].sort_values("d_km")
        ax.plot(sub["d_km"], sub[ycol], 'o-', color=color, ms=4, label=f"{case} turbulence")
        if ci and ycol == "R_mc":
            ax.fill_between(sub["d_km"], sub["R_mc_ci95_lo"], sub["R_mc_ci95_hi"], color=color, alpha=0.15)
    if det_reference is not None:
        ax.plot(det_reference["d_km"], det_reference["R"], 'k--', lw=1.3, label="Deterministic FSO")
    if logy:
        ax.set_yscale('log')
    ax.set_xlabel("One-way FSO distance $d$ (km)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, name)
    plt.show()
    return fig


def plot_channel_loss_with_fading(distances_main=DISTANCES_MAIN, seed=1):
    rng_loss = np.random.default_rng(seed)
    fig, ax = plt.subplots(figsize=(7, 5))
    det_loss = channel_loss_dB(distances_main)
    ax.plot(distances_main, det_loss, 'k--', lw=1.3, label="Deterministic loss (Model A)")
    for case, cn2 in TURBULENCE_CASES.items():
        lo5, hi95, med = [], [], []
        for d in distances_main:
            s2 = rytov_variance(d, cn2)
            ag, bg = gamma_gamma_params(s2)
            h = sample_gamma_gamma(ag, bg, 20_000, rng_loss)
            fading_loss_db = -10 * np.log10(h)
            lo5.append(np.percentile(fading_loss_db, 5))
            hi95.append(np.percentile(fading_loss_db, 95))
            med.append(np.percentile(fading_loss_db, 50))
        total_med = det_loss + np.array(med)
        total_lo = det_loss + np.array(lo5)
        total_hi = det_loss + np.array(hi95)
        ax.plot(distances_main, total_med, '-', color=CASE_COLORS[case], label=f"{case}: median total loss")
        ax.fill_between(distances_main, total_lo, total_hi, color=CASE_COLORS[case], alpha=0.15)
    ax.set_xlabel("One-way FSO distance $d$ (km)")
    ax.set_ylabel("Channel loss (dB)")
    ax.set_title("Channel loss: deterministic + Gamma-Gamma fading (5th-95th pct band)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "gg_channel_loss_vs_distance")
    plt.show()
    return fig


def plot_all_gamma_gamma_figures(df_det, df_gg):
    plot_gg_vs_distance(df_gg, "R_mc", "Ensemble key rate $\\langle R\\rangle$", "Gamma-Gamma: SKR vs distance (Monte Carlo)",
                         "gg_skr_vs_distance", det_reference=df_det, ci=True)
    plot_gg_vs_distance(df_gg, "Qxx", "$\\langle Q_{xx}\\rangle$", "Gamma-Gamma: $Q_{xx}$ vs distance", "gg_qxx_vs_distance", logy=False)
    plot_gg_vs_distance(df_gg, "Exx", "$\\langle E_{xx}\\rangle$", "Gamma-Gamma: $E_{xx}$ vs distance", "gg_exx_vs_distance", logy=False)
    plot_gg_vs_distance(df_gg, "s1", "$\\langle s_1\\rangle$", "Gamma-Gamma: $s_1$ vs distance", "gg_s1_vs_distance", logy=False)
    plot_channel_loss_with_fading()


# ---------------------------------------------------------------------------
# Part 9: optimisation-curve figures
# ---------------------------------------------------------------------------
def plot_optimal_param_vs_distance(df_det, df_gg, param, ylabel, title, name):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(df_det["d_km"], df_det[param], 'k--', lw=1.3, label="Deterministic FSO")
    for case, color in CASE_COLORS.items():
        sub = df_gg[df_gg["turbulence_case"] == case].sort_values("d_km")
        ax.plot(sub["d_km"], sub[param], 'o-', color=color, ms=4, label=f"{case} (ensemble-optimal)")
    ax.set_xlabel("One-way FSO distance $d$ (km)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, name)
    plt.show()
    return fig


def plot_epsilon_improvement(df_gg, distances_main=DISTANCES_MAIN, e_a=E_A,
                              eps_naive_fail=0.05, eps_naive_mild=0.01):
    R_opt_by_case, R_fail_by_case, R_mild_by_case = {}, {}, {}
    for case in TURBULENCE_CASES:
        sub = df_gg[df_gg["turbulence_case"] == case].sort_values("d_km")
        R_opt, R_fail, R_mild = [], [], []
        for _, row in sub.iterrows():
            s2 = rytov_variance(row["d_km"], TURBULENCE_CASES[case])
            R_opt.append(row["R_quad"])
            R_fail.append(gauss_laguerre_ensemble_rate(eps_naive_fail, row["mu"], row["d_km"], e_a, s2, n_points=40))
            R_mild.append(gauss_laguerre_ensemble_rate(eps_naive_mild, row["mu"], row["d_km"], e_a, s2, n_points=40))
        R_opt_by_case[case], R_fail_by_case[case], R_mild_by_case[case] = R_opt, R_fail, R_mild

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for case, color in CASE_COLORS.items():
        axes[0].plot(distances_main, R_opt_by_case[case], 'o-', color=color, ms=4, label=f"{case}: optimal $\\epsilon^*$")
        axes[0].plot(distances_main, np.clip(R_fail_by_case[case], 1e-30, None), 'o--', color=color, ms=3, alpha=0.5,
                     label=f"{case}: naive $\\epsilon=0.05$ (fails)")
    axes[0].set_yscale('log')
    axes[0].set_xlabel("One-way FSO distance $d$ (km)")
    axes[0].set_ylabel("Ensemble key rate $\\langle R\\rangle$")
    axes[0].set_title("Optimal $\\epsilon^*$ vs a poorly-chosen $\\epsilon=0.05$")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.3)
    for case, color in CASE_COLORS.items():
        ratio = np.array(R_opt_by_case[case]) / np.array(R_mild_by_case[case])
        axes[1].plot(distances_main, ratio, 'o-', color=color, ms=4, label=case)
    axes[1].axhline(1.0, color='gray', ls=':', lw=1)
    axes[1].set_xlabel("One-way FSO distance $d$ (km)")
    axes[1].set_ylabel("Improvement factor $R(\\epsilon^*)/R(\\epsilon=0.01)$")
    axes[1].set_title("SKR improvement from optimising $\\epsilon$\n(mild baseline $\\epsilon=0.01$)")
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "epsilon_skr_improvement_factor")
    plt.show()
    return fig, (R_opt_by_case, R_fail_by_case, R_mild_by_case)


def plot_instantaneous_map(df_inst, value_col, ylabel, title, name):
    pivot = df_inst.pivot(index="d_km", columns="percentile", values=value_col)
    fig, ax = plt.subplots(figsize=(7, 5.5))
    im = ax.imshow(pivot.values, aspect='auto', cmap='viridis', origin='lower')
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([f"{c:g}" for c in pivot.columns])
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels([f"{d:g}" for d in pivot.index], fontsize=8)
    ax.set_xlabel("Turbulence fading percentile")
    ax.set_ylabel("One-way FSO distance $d$ (km)")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=ylabel)
    fig.tight_layout()
    save_fig(fig, name)
    plt.show()
    return fig


def plot_joint_optimization_landscape(d_demo=6.0, cn2=None, e_a=E_A, n_grid=50, n_points=20):
    from src.fso.config import CN2_MODERATE
    if cn2 is None:
        cn2 = CN2_MODERATE
    s2_demo = rytov_variance(d_demo, cn2)
    eps_grid = np.linspace(*EPS_BOUNDS, n_grid)
    mu_grid = np.linspace(*MU_BOUNDS, n_grid)
    R_land = np.zeros((n_grid, n_grid))
    for i, e in enumerate(eps_grid):
        for j, m in enumerate(mu_grid):
            R_land[i, j] = gauss_laguerre_ensemble_rate(e, m, d_demo, e_a, s2_demo, n_points=n_points)
    opt_demo = joint_optimize_ensemble(d_demo, e_a, s2_demo, seed=0)
    print(opt_demo)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    vmax = max(R_land.max(), opt_demo["R"])
    im = ax.pcolormesh(mu_grid, eps_grid, R_land, shading='auto', cmap='viridis', vmin=0.0, vmax=vmax)
    ax.plot(opt_demo["mu"], opt_demo["eps"], 'r*', ms=18, mec='white', mew=1.2, label='joint optimum')
    ax.set_xlabel("Signal intensity $\\mu'$")
    ax.set_ylabel("Sending probability $\\epsilon$")
    ax.set_ylim(0, 0.3)
    ax.set_title(f"Joint optimisation landscape (d={d_demo} km)\n(negative-rate region clipped to 0 for visibility; $\\epsilon$ axis zoomed)")
    fig.colorbar(im, ax=ax, label="Ensemble key rate $\\langle R\\rangle$ (clipped at 0)")
    ax.legend(fontsize=9, loc='upper right')
    fig.tight_layout()
    save_fig(fig, f"joint_optimization_heatmap_d{d_demo}km")
    plt.show()
    return fig, opt_demo


def plot_quadrature_validation(df_quad_val):
    fig, ax = plt.subplots(figsize=(7, 5))
    markers = {20: 'o', 40: 's', 80: '^'}
    for n in (20, 40, 80):
        sub = df_quad_val[df_quad_val.n_points == n]
        ax.plot(sub.d_km, sub.rel_err * 100, marker=markers[n], label=f"{n}-point GL")
    ax.axhline(1.0, color='gray', ls=':', label='1% target')
    ax.set_yscale('log')
    ax.set_xlabel("One-way FSO distance $d$ (km)")
    ax.set_ylabel("Relative error vs Monte Carlo (%)")
    ax.set_title("Gauss-Laguerre quadrature validation")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "quadrature_validation_error")
    plt.show()
    return fig


def plot_monte_carlo_convergence(n_arr, running_mean, true_val, d_demo, case_label="moderate turbulence"):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(n_arr, running_mean, color='#1d3557')
    ax.axhline(true_val, color='#e76f51', ls='--', label='Gauss-Laguerre reference')
    ax.set_xscale('log')
    ax.set_xlabel("Number of Monte Carlo samples $N_{MC}$")
    ax.set_ylabel("Running-mean ensemble key rate")
    ax.set_title(f"Monte Carlo convergence (d={d_demo} km, {case_label})")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "monte_carlo_convergence")
    plt.show()
    return fig
