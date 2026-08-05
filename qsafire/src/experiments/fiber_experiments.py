"""Part I experiments: illustrations, validation tables, and Figure 1/2
reproductions (notebook Sections 3, 4, 5, 6, 7, 8, 9, 10, 11, 12).
"""

import numpy as np
import pandas as pd

from src.engine import (
    ETA_D, P_DARK,
    gain_X, error_X, S_Z_gain, E_Z_qber,
    s1_decoy_lower_bound, e1ph_upper_bound,
)
from src.simulation.fiber_channel import (
    ALPHA, channel_eta, true_single_photon_yield, sns_key_rate, optimize_key_rate,
)
from src.visualization.fiber_plots import (
    plot_gain_and_error_vs_intensity, plot_z_qber_vs_epsilon, plot_figure1, plot_figure2,
)


def demo_channel_scaling():
    """Section 3: quick sanity check / illustration of the sqrt-scaling
    advantage."""
    L_demo = np.array([0, 100, 200, 400, 600, 800])
    eta_tfqkd = channel_eta(L_demo)
    eta_pp = ETA_D * 10 ** (-ALPHA * L_demo / 10)  # ordinary point-to-point transmittance
    df = pd.DataFrame({
        'L (km)': L_demo,
        'eta_TFQKD(L) [~sqrt scaling]': eta_tfqkd,
        'eta_point-to-point(L) [linear scaling]': eta_pp,
    })
    return df


def demo_gain_and_error(show_plot=True):
    """Section 4: gain/error illustration at L=300 km."""
    if show_plot:
        plot_gain_and_error_vs_intensity(300)


def demo_z_qber(show_plot=True):
    """Section 5: Z-basis QBER vs eps illustration at L=300 km, mu'=0.3."""
    if show_plot:
        plot_z_qber_vs_epsilon(300, 0.3)


def validate_s1_decoy_estimate():
    """Section 6: compare the decoy-estimated s1 against the true (model)
    single-photon yield Y1(L), as L is varied."""
    rows = []
    for L in [0, 100, 200, 300, 400, 500, 600, 700, 800]:
        eta = channel_eta(L)
        mu2 = 0.30
        mu1 = 0.06
        s0 = 2.0 * P_DARK
        S1 = gain_X(mu1, eta)
        S2 = gain_X(mu2, eta)
        s1_est = s1_decoy_lower_bound(mu1, mu2, S1, S2, s0)
        s1_true = true_single_photon_yield(eta)
        rows.append((L, eta, s1_true, s1_est, abs(s1_true - s1_est) / s1_true))

    df_s1 = pd.DataFrame(rows, columns=["L (km)", "eta(L)", "s1_true", "s1_decoy_estimate", "relative_error"])
    return df_s1


def demo_e1ph_recovery():
    """Section 7: recovered e1^ph vs the "true" injected misalignment e_a, at
    L=400 km."""
    eta_demo = channel_eta(400)
    mu1, mu2 = 0.06, 0.30
    s0 = 2.0 * P_DARK
    S1 = gain_X(mu1, eta_demo)
    S2 = gain_X(mu2, eta_demo)
    s1_est = s1_decoy_lower_bound(mu1, mu2, S1, S2, s0)

    ea_values = np.linspace(0.0, 0.45, 10)
    recovered = []
    for ea in ea_values:
        E1 = error_X(mu1, eta_demo, ea)
        recovered.append(e1ph_upper_bound(mu1, S1, E1, s0, s1_est))

    return pd.DataFrame({"injected e_a": ea_values, "recovered e1^ph (Eq. 45)": recovered})


def quick_check_key_rate():
    """Section 8 quick check."""
    R = sns_key_rate(0.05, 0.3, 200, 0.15)
    print("R(L=200km, ea=0.15, eps=0.05, mu'=0.3) =", R)
    return R


def quick_check_optimizer():
    """Section 9 quick check."""
    R_opt, eps_opt, mu_opt = optimize_key_rate(300, 0.25)
    print(f"Optimised @ L=300km, e_a=0.25:  R={R_opt:.4e}, eps={eps_opt:.4e}, mu'={mu_opt:.4e}")
    return R_opt, eps_opt, mu_opt


def figure1_sweep(show_plot=True):
    """Section 10: reproduce Figure 1 (key rate vs distance for several
    e_a)."""
    distances = np.linspace(1, 900, 46)  # km
    ea_values_fig1 = [0.00, 0.15, 0.25, 0.35, 0.45]

    fig1_results = {}
    for ea in ea_values_fig1:
        Rs, eps_list, mu_list = [], [], []
        for L in distances:
            R, e, m = optimize_key_rate(L, ea)
            Rs.append(R)
            eps_list.append(e)
            mu_list.append(m)
        fig1_results[ea] = {
            "R": np.array(Rs),
            "eps": np.array(eps_list),
            "mu": np.array(mu_list),
        }

    print("Optimisation sweep for Figure 1 complete.")

    if show_plot:
        plot_figure1(distances, fig1_results, ea_values_fig1)

    print("Approximate maximum secure distance for each misalignment level:")
    for ea in ea_values_fig1:
        R = fig1_results[ea]["R"]
        secure_L = distances[R > 0]
        max_L = secure_L.max() if len(secure_L) > 0 else 0.0
        print(f"  e_a = {ea*100:4.0f}%  ->  max secure distance ~ {max_L:.0f} km")

    return distances, fig1_results, ea_values_fig1


def figure2_sweep(L_fixed=500.0, show_plot=True):
    """Section 11: reproduce Figure 2 (key rate vs misalignment at fixed
    L=500km)."""
    ea_sweep = np.linspace(0.0, 0.35, 36)

    fig2_R, fig2_eps, fig2_mu = [], [], []
    for ea in ea_sweep:
        R, e, m = optimize_key_rate(L_fixed, ea)
        fig2_R.append(R)
        fig2_eps.append(e)
        fig2_mu.append(m)

    fig2_R = np.array(fig2_R)

    if show_plot:
        plot_figure2(ea_sweep, fig2_R, L_fixed)

    secure_ea = ea_sweep[fig2_R > 0]
    if len(secure_ea):
        print(f"Largest tolerable misalignment error at L={L_fixed:.0f} km with positive key rate: "
              f"~{secure_ea.max()*100:.1f}%")
    else:
        print("No positive key rate found.")

    return ea_sweep, fig2_R, fig2_eps, fig2_mu


def full_diagnostics(L, e_a, decoy_ratio=0.2):
    """Section 12: run the optimiser and then recompute every intermediate
    quantity for reporting."""
    R_opt, eps_opt, mu_opt = optimize_key_rate(L, e_a, decoy_ratio)

    eta = channel_eta(L)
    mu2 = mu_opt
    mu1 = decoy_ratio * mu2
    s0 = 2.0 * P_DARK

    S1 = gain_X(mu1, eta)
    S2 = gain_X(mu2, eta)
    E1 = error_X(mu1, eta, e_a)

    s1 = max(s1_decoy_lower_bound(mu1, mu2, S1, S2, s0), 1e-15)
    e1ph = float(np.clip(e1ph_upper_bound(mu1, S1, E1, s0, s1), 0.0, 0.5))

    Sz = S_Z_gain(eps_opt, mu_opt, eta)
    Ez = float(np.clip(E_Z_qber(eps_opt, mu_opt, eta), 1e-12, 0.5))

    return {
        "L (km)": L, "e_a": e_a,
        "eps*": eps_opt, "mu'*": mu_opt,
        "s1": s1, "e1^ph": e1ph,
        "S_Z": Sz, "E^Z": Ez,
        "R (bits/pulse)": R_opt,
    }


def validation_tables():
    """Section 12: full validation table plus the focused e_a=0.25 view."""
    rows = []
    for L in [100, 300, 500, 700]:
        for ea in [0.0, 0.15, 0.25, 0.35]:
            rows.append(full_diagnostics(L, ea))

    df_validation = pd.DataFrame(rows)
    pd.set_option('display.float_format', lambda v: f'{v:.6g}')

    df_focus = df_validation[df_validation["e_a"] == 0.25].reset_index(drop=True)
    print("Representative optimised parameters at e_a = 25% misalignment:\n")
    print(df_focus.to_string(index=False))

    return df_validation, df_focus
