# snstfqkd/optimization/optimizer.py
#
# Verbatim move from the original top-level optimizer.py. The only
# structural change: the duplicated _get_S_Z() helper (identical in
# the original key_rate.py and optimizer.py) now lives once in
# security/bounds.py and is imported here — no formula changed.
#
# PATCH NOTES (carried over from original):
#   _evaluate() and optimize_all() use z_window_yield for S_Z.
#   Physical justification: same as key_rate.py patch.
#   Reference: Wang 2018, Eq. (4), Step 3 Z-window definition.
import numpy as np
from dataclasses import dataclass
from typing import Tuple

from scipy.optimize import differential_evolution

from snstfqkd.security.decoy_analysis import DecoyStateAnalyzer
from snstfqkd.security.key_rate import SNSKeyRateCalculator, binary_entropy
from snstfqkd.security.bounds import get_S_Z as _get_S_Z


@dataclass
class OptimizationResult:
    key_rate: float
    epsilon: float
    mu_signal: float
    mu1: float
    mu2: float
    s1: float
    e1_ph: float
    method: str = "unknown"


class SNSOptimizer:
    """
    Jointly optimises epsilon, mu', mu1, mu2 to maximise the asymptotic
    SNS key rate.

    Uses z_window_yield() for S_Z to correctly model the Z-window
    counting rate (Wang 2018 Eq. 4).
    """

    def __init__(
        self,
        channel,
        f: float = 1.1,
        epsilon_bounds: Tuple[float, float] = (1e-3, 0.5),
        mu_signal_bounds: Tuple[float, float] = (1e-3, 0.5),
        mu1_bounds: Tuple[float, float] = (1e-4, 0.3),
        mu2_bounds: Tuple[float, float] = (0.05, 0.5),
        de_maxiter: int = 300,
    ):
        self.channel = channel
        self.calculator = SNSKeyRateCalculator(f=f)
        self.analyzer = DecoyStateAnalyzer()
        self.epsilon_bounds = epsilon_bounds
        self.mu_signal_bounds = mu_signal_bounds
        self.mu1_bounds = mu1_bounds
        self.mu2_bounds = mu2_bounds
        self.de_maxiter = de_maxiter

    def _evaluate(self, distance, epsilon, mu_signal, mu1, mu2) -> float:
        """
        Evaluate key rate for one parameter set.

        S_Z uses z_window_yield (one-arm BS) not intensity_yield
        (two-arm interference). This is the correct observable for Wang
        Eq. (4).
        """
        if mu1 >= mu2:
            return 0.0
        ch = self.channel
        s0    = ch.vacuum_yield()
        S_mu1 = ch.intensity_yield(mu1, distance)    # X-window
        S_mu2 = ch.intensity_yield(mu2, distance)    # X-window
        E_mu1 = ch.qber_x_basis(mu1, distance)      # X-window QBER
        E_Z   = ch.qber_z_basis()
        S_Z   = _get_S_Z(ch, mu_signal, distance)
        decoy = self.analyzer.three_intensity_decoy(mu1, mu2, s0, S_mu1, S_mu2, E_mu1)
        if not decoy.valid:
            return 0.0
        return self.calculator.key_rate(epsilon, mu_signal, decoy.s1, decoy.e1_ph, S_Z, E_Z)

    def optimize_epsilon(self, distance, mu_signal, mu1, mu2, steps=50) -> float:
        best_eps, best_rate = self.epsilon_bounds[0], 0.0
        for eps in np.linspace(*self.epsilon_bounds, steps):
            r = self._evaluate(distance, eps, mu_signal, mu1, mu2)
            if r > best_rate:
                best_rate, best_eps = r, eps
        return best_eps

    def optimize_mu_signal(self, distance, epsilon, mu1, mu2, steps=50) -> float:
        best_mu, best_rate = self.mu_signal_bounds[0], 0.0
        for mu in np.linspace(*self.mu_signal_bounds, steps):
            r = self._evaluate(distance, epsilon, mu, mu1, mu2)
            if r > best_rate:
                best_rate, best_mu = r, mu
        return best_mu

    def optimize_decoy_parameters(self, distance, epsilon, mu_signal, steps=30) -> Tuple[float, float]:
        best_mu1, best_mu2, best_rate = self.mu1_bounds[0], self.mu2_bounds[0], 0.0
        for mu1 in np.linspace(*self.mu1_bounds, steps):
            for mu2 in np.linspace(*self.mu2_bounds, steps):
                if mu1 >= mu2:
                    continue
                r = self._evaluate(distance, epsilon, mu_signal, mu1, mu2)
                if r > best_rate:
                    best_rate, best_mu1, best_mu2 = r, mu1, mu2
        return best_mu1, best_mu2

    def optimize_all(
        self,
        distance: float,
        method: str = "differential_evolution",
        grid_steps: int = 20,
        de_maxiter: int = None,
        de_popsize: int = 12,
        seed: int = 42,
    ) -> OptimizationResult:
        """Jointly optimise all four free parameters.

        Uses corrected _evaluate() which calls z_window_yield for S_Z.
        """
        if de_maxiter is None:
            de_maxiter = self.de_maxiter

        bounds = [
            self.epsilon_bounds,
            self.mu_signal_bounds,
            self.mu1_bounds,
            self.mu2_bounds,
        ]

        if method == "differential_evolution":
            def neg_rate(x):
                return -self._evaluate(distance, *x)

            res = differential_evolution(
                neg_rate, bounds,
                maxiter=de_maxiter, popsize=de_popsize,
                seed=seed, tol=1e-12,
                mutation=(0.5, 1.5), recombination=0.7,
            )
            eps, mu_s, mu1, mu2 = res.x
            rate = -res.fun

        else:  # grid
            eps_g  = np.linspace(*self.epsilon_bounds, grid_steps)
            mu_s_g = np.linspace(*self.mu_signal_bounds, grid_steps)
            mu1_g  = np.linspace(*self.mu1_bounds, grid_steps // 2)
            mu2_g  = np.linspace(*self.mu2_bounds, grid_steps // 2)
            eps, mu_s, mu1, mu2, rate = (
                self.epsilon_bounds[0], self.mu_signal_bounds[0],
                self.mu1_bounds[0], self.mu2_bounds[0], 0.0,
            )
            for e in eps_g:
                for ms in mu_s_g:
                    for m1 in mu1_g:
                        for m2 in mu2_g:
                            if m1 >= m2:
                                continue
                            r = self._evaluate(distance, e, ms, m1, m2)
                            if r > rate:
                                rate, eps, mu_s, mu1, mu2 = r, e, ms, m1, m2

        # Recover decoy quantities for the winning parameter set
        s0    = self.channel.vacuum_yield()
        S_mu1 = self.channel.intensity_yield(mu1, distance)
        S_mu2 = self.channel.intensity_yield(mu2, distance)
        E_mu1 = self.channel.qber_x_basis(mu1, distance)
        decoy = self.analyzer.three_intensity_decoy(mu1, mu2, s0, S_mu1, S_mu2, E_mu1)

        return OptimizationResult(
            key_rate=max(0.0, rate),
            epsilon=float(eps), mu_signal=float(mu_s),
            mu1=float(mu1), mu2=float(mu2),
            s1=decoy.s1, e1_ph=decoy.e1_ph,
            method=method,
        )
