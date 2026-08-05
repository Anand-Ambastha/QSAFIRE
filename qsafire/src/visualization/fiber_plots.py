"""Plotting routines for Part I (fibre SNS-TF-QKD reproduction).

Every figure below preserves the exact colors, labels, legends, scales,
annotations, titles and layout of the corresponding notebook cell.
"""

import numpy as np
import matplotlib.pyplot as plt

from src.simulation.fiber_channel import channel_eta
from src.engine import gain_X, error_X, E_Z_qber


def plot_channel_scaling_illustration():
    """Section 3 quick sanity-check table is numeric only (see
    fiber_experiments.demo_channel_scaling); no plot in that cell."""
    pass


def plot_gain_and_error_vs_intensity(L_demo_km=300):
    """Section 4 illustration: gain and error vs decoy intensity at a fixed
    distance."""
    eta_demo = channel_eta(L_demo_km)
    mus = np.linspace(0.001, 0.6, 100)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(mus, gain_X(mus, eta_demo))
    axes[0].set_xlabel("decoy intensity $\\mu$")
    axes[0].set_ylabel("$S_\\mu$ (gain)")
    axes[0].set_title(f"Observed gain vs decoy intensity, L={L_demo_km} km")
    for ea in [0.0, 0.15, 0.25, 0.35, 0.45]:
        axes[1].plot(mus, error_X(mus, eta_demo, ea), label=f"$e_a$={ea}")
    axes[1].set_xlabel("decoy intensity $\\mu$")
    axes[1].set_ylabel("$E_\\mu^X$ (error rate)")
    axes[1].set_title(f"X-basis error vs decoy intensity, L={L_demo_km} km")
    axes[1].legend(fontsize=8)
    plt.tight_layout()
    plt.show()
    return fig


def plot_z_qber_vs_epsilon(L_demo_km=300, mu_prime=0.3):
    """Section 5 illustration: E^Z stays extremely small regardless of
    misalignment e_a, because misalignment never enters this formula."""
    eta_demo = channel_eta(L_demo_km)
    eps_demo = np.linspace(0.001, 0.3, 100)
    Ez_vals = [E_Z_qber(e, mu_prime, eta_demo) for e in eps_demo]
    fig = plt.figure(figsize=(5.5, 4))
    plt.plot(eps_demo, Ez_vals)
    plt.xlabel("sending probability $\\epsilon$")
    plt.ylabel("$E^Z$")
    plt.title(f"Z-basis QBER vs $\\epsilon$ (L={L_demo_km} km, $\\mu'$={mu_prime})\n"
              "-> independent of single-photon-interference misalignment $e_a$")
    plt.tight_layout()
    plt.show()
    return fig


def plot_figure1(distances, fig1_results, ea_values_fig1):
    """Section 10: Reproduction of Figure 1 (key rate vs. distance)."""
    fig = plt.figure(figsize=(7.5, 5.5))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(ea_values_fig1)))
    for ea, c in zip(ea_values_fig1, colors):
        R = fig1_results[ea]["R"]
        mask = R > 0
        plt.plot(distances[mask], np.log10(R[mask]), label=f"$e_a$ = {ea*100:.0f}%", color=c, lw=2)

    plt.xlabel("Distance between Alice and Bob (km)")
    plt.ylabel(r"$\log_{10}$(key rate per pulse)")
    plt.title("Reproduction of Fig. 1: SNS-TFQKD key rate vs. distance\n"
              "(operational protocol + Protocol-4 reduction + decoy analysis, Eqs. 4, 44, 45)")
    plt.ylim(-13, 0)
    plt.xlim(0, 900)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()
    return fig


def plot_figure2(ea_sweep, fig2_R, L_fixed=500.0):
    """Section 11: Reproduction of Figure 2 (key rate vs. misalignment @
    500 km)."""
    fig = plt.figure(figsize=(6.5, 5))
    mask = fig2_R > 0
    plt.plot(ea_sweep[mask], np.log10(fig2_R[mask]), 'b-o', markersize=3)
    plt.xlabel("Misalignment error $e_a$")
    plt.ylabel(r"$\log_{10}$(key rate per pulse)")
    plt.title(f"Reproduction of Fig. 2: key rate vs. misalignment\n(Alice-Bob distance fixed at {L_fixed:.0f} km)")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
    return fig
