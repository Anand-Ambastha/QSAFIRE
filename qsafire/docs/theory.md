# Theory

Physics and mathematics behind each module, in the order the link budget is
built up: protocol core -> channel models -> satellite link budget layers ->
optimization -> multi-link generalization. Implementation details (function
signatures, parameters) are in `docs/api.md`.

---

## 1. The SNS-TF-QKD protocol (`src/engine.py`)

Sending-or-Not-Sending Twin-Field QKD (Wang, Yu & Hu, Phys. Rev. A **98**,
062323, 2018) is a twin-field protocol: Alice and Bob each send weak
coherent pulses toward a central untrusted relay, Charlie, who performs an
interference measurement and announces successful clicks. Its headline
advantage over ordinary repeaterless QKD is a **square-root scaling** of the
key rate with total channel transmittance, rather than the linear scaling of
older protocols — because each arm contributes to Charlie's click
probability additively (through the combined optical field), not
multiplicatively.

**Basis structure.** In each round, Alice and Bob each independently decide,
with probability `eps`, to send a phase-randomized coherent pulse (the
"send" decision) in the **Z window**, used for the raw key; with the
complementary probability they instead send a *test* pulse in the **X
window**, used to estimate the channel's error rate and (via decoy states)
the single-photon yield and phase-error rate.

**Photon-number weights.** For a coherent state of mean photon number `mu`,
the probability of `n` photons follows the Poisson distribution: `p0(mu) =
e^{-mu}`, `p1(mu) = mu*e^{-mu}`, `p2(mu) = (mu^2/2)*e^{-mu}`.

**Gain and QBER.** The X-window gain `S_mu` (probability Charlie reports a
successful detection given intensity `mu` and channel transmittance `eta`)
and its associated QBER `E_mu` follow directly from the detector click
model: a real click occurs with probability `1-e^{-mu*eta}`, a false
(dark-count) click occurs with probability proportional to `p_dark`, and the
error rate blends the two through the reported optical misalignment error
`e_a`.

**Z-window gain/QBER** (`S_Z_gain`, `E_Z_qber`) follow the same click model
at the "send" intensity `mu'`, but weighted by the sending probability
`eps` — and, per the protocol's own well-known feature, the Z-basis QBER is
dominated only by dark counts (not by optical misalignment), since Z-basis
bits are only kept when exactly one party sent.

**Decoy-state estimation.** Because Alice and Bob cannot directly measure
the single-photon yield `s1` or single-photon phase-error rate `e1ph` (they
only observe pulse-averaged gains/QBERs across the Poisson-mixed photon
numbers actually sent), the protocol uses two decoy intensities (`mu1 <
mu2`) and linear-programming-style bounds:

- `s1_decoy_lower_bound(mu1, mu2, S_mu1, S_mu2, s0)` — a lower bound on the
  true single-photon yield, from the two measured gains and the (separately
  bounded) vacuum yield `s0`.
- `e1ph_upper_bound(mu1, S_mu1, E_mu1, s0, s1)` — an upper bound on the
  single-photon phase-error rate, from the measured QBER at the weaker decoy
  intensity.

These are the standard two-decoy-state bounds used throughout the
decoy-state QKD literature, applied here to the SNS protocol's own gain/QBER
definitions.

**The asymptotic key rate.**

```
R = 2*eps*(1-eps)*mu'*e^{-mu'}*s1*(1 - h2(e1ph)) - Sz*f_ec*h2(Ez)
```

where `h2` is the binary Shannon entropy. The first term is the raw
single-photon key generation rate (Z-window sifted fraction, weighted by the
single-photon yield and privacy-amplification cost `1-h2(e1ph)`); the second
term is the error-correction cost (`f_ec` is the inefficiency factor above
the Shannon limit a real error-correction code needs). `R` can be negative
when the channel is too lossy or too noisy for the estimated bounds to
support a positive key — this is a real, physically meaningful "no secure
key possible" result, not a numerical artifact, and every downstream module
in this package clips negative `R` to zero only when *accumulating* key
volume (a negative rate cannot be "spent back"), while preserving the raw
signed value everywhere it is used for optimization or diagnosis.

---

## 2. Fiber channel (`src/simulation/fiber_channel.py`)

The simplest channel model: transmittance decays exponentially with
distance per the standard fiber attenuation coefficient (`alpha`, dB/km),
scaled by a fixed detector efficiency `eta_d`:

```
eta(L) = eta_d * 10^(-alpha*L/10)
```

`sns_key_rate` and `optimize_key_rate` are thin wrappers around
`src.engine.sns_key_rate_eta`/its optimizer, parameterized by distance
instead of transmittance directly.

---

## 3. FSO channel with Gamma-Gamma turbulence (`src/fso/`)

**Deterministic loss.** Free-space paths lose optical power through beam
spreading and atmospheric extinction; `fso_deterministic_eta` models this
with the same exponential-attenuation-coefficient form as fiber, but with an
FSO-representative attenuation coefficient.

**Turbulence and the Rytov variance.** Atmospheric turbulence causes random
refractive-index fluctuations along the propagation path, which in turn
cause the received irradiance to fluctuate ("scintillation"). The strength
of this effect for weak-to-moderate turbulence is captured by the **Rytov
variance**,

```
sigma_R^2 = 1.23 * Cn2 * k^(7/6) * L^(11/6)     (plane wave, horizontal path)
```

where `Cn2` is the refractive-index structure parameter (m^(-2/3)), `k =
2*pi/lambda` is the optical wavenumber, and `L` is the path length. Larger
`sigma_R^2` means stronger turbulence.

**Gamma-Gamma fading.** The Gamma-Gamma distribution is the standard model
for irradiance fluctuations spanning weak-to-strong turbulence, derived as
the product of two independent Gamma-distributed random variables
representing large-scale and small-scale eddies. Its two shape parameters
`(alpha_g, beta_g)` are related to the Rytov variance through the standard
closed-form expressions (`gamma_gamma_params`), and turbulence is classified
as weak/moderate/strong (`classify_turbulence_regime`) based on where
`sigma_R^2` falls relative to the standard weak-fluctuation validity bound
(`validity_flag`).

**Ensemble averaging.** Because the key rate formula is nonlinear in `eta`,
the *mean* key rate under fading is **not** the key rate evaluated at the
mean transmittance — it must be averaged over the Gamma-Gamma fading
distribution. This package computes that average two independent ways,
cross-validated against each other:

- **Monte Carlo** (`monte_carlo.py`): draw many samples `h` from the
  Gamma-Gamma distribution, evaluate the (fast, closed-form) key rate at
  `eta = eta_det * h` for each sample, and average.
- **Gauss-Laguerre quadrature** (`quadrature.py`): since the Gamma-Gamma PDF
  is a product of two Gamma densities, the ensemble average can be computed
  as a small number of weighted quadrature-node evaluations instead of many
  random samples — much faster for repeated evaluation inside an optimizer,
  at the cost of a small (and directly measured) quadrature-truncation
  error.

**Optimization** (`optimization.py`) triangulates local grid search, global
Differential Evolution, and Basin Hopping to jointly optimize `(eps, mu')`
under ensemble averaging, since the ensemble-averaged rate surface can have
local optima that a purely local optimizer would miss.

---

## 4. Ground-to-satellite geometry (`src/satellite/geometry.py`)

**Geodesy.** Ground-station positions are given in WGS-84 geodetic
coordinates (latitude, longitude, altitude) and converted to Earth-Centered
Earth-Fixed (ECEF) Cartesian coordinates using the standard WGS-84
ellipsoid-of-revolution formula (`geodetic_to_ecef`), accounting for the
Earth's flattening.

**Sun-synchronous orbit design.** A sun-synchronous orbit is one whose
orbital plane precesses (due to the Earth's oblateness, characterized by the
`J2` zonal harmonic) at exactly the rate the Earth orbits the Sun, so the
satellite crosses a given latitude at the same local solar time on every
pass. Rather than assume a textbook inclination, `sun_sync_orbit` **solves**
for the inclination that satisfies the secular nodal-precession condition

```
d(RAAN)/dt = -1.5 * n0 * J2 * (Re/p)^2 * cos(i) = omega_sun
```

at the requested altitude, together with the corresponding first-order J2
secular rates for the argument of latitude — the standard perturbation
result for a near-circular low-Earth orbit.

**Propagation.** The satellite's position is propagated in the Earth-Centered
Inertial (ECI) frame from its orbital elements (`eci_position`), then
rotated into ECEF using Greenwich Mean Sidereal Time (GMST, computed from
the Julian date via the standard IAU low-precision series in `gmst_rad`) to
account for Earth's rotation (`eci_to_ecef`).

**Look angles.** Given the satellite's ECEF position and a ground station's
ECEF position, `look_angles` projects the relative vector into the site's
local East-North-Up (ENU) frame to get elevation, azimuth, and slant range —
the quantities that drive every downstream loss calculation.

**Pass geometry.** A specific dual-visibility pass (both ground stations
above a minimum elevation simultaneously) is engineered by choosing the
orbital plane's right ascension of the ascending node (RAAN) and initial
argument of latitude so the satellite crosses the two stations'
great-circle midpoint at a chosen epoch, then propagating forward/backward
in time to find the common-visibility window.

---

## 5. Static atmospheric transmittance (`src/atmosphere/transmittance.py`)

Atmospheric extinction is modeled by Beer's law, `eta_atm = exp(-tau)`,
where the total **zenith optical depth** `tau` is the sum of three
independent extinction mechanisms:

- **Rayleigh (molecular) scattering** (`tau_rayleigh_zenith`): scattering
  off air molecules, following the standard wavelength-dependent empirical
  formula, scaled linearly with surface pressure (thinner air at higher
  station altitude means less Rayleigh extinction).
- **Aerosol (Mie) extinction** (`tau_aerosol_zenith`): scattering and
  absorption by suspended particulates, parameterized by the column Aerosol
  Optical Depth (AOD) at a reference wavelength and scaled to the operating
  wavelength via the Angstrom power law, `tau(lambda) = AOD(lambda0) *
  (lambda/lambda0)^(-alpha)`, where `alpha` (the Angstrom exponent)
  characterizes the aerosol particle-size distribution. A separate
  visibility-based cross-check (`tau_aerosol_visibility`, via the Kim/Kruse
  size-distribution parameter) provides an independent estimate from a
  different observable (meteorological visibility) for validation.
- **Gaseous absorption**: a small fixed residual for absorption bands not
  captured by the Rayleigh/aerosol terms.

**Slant-path scaling.** The zenith optical depth applies to a path straight
up; for a path at elevation angle `el` (zenith angle `90-el`), the optical
path length through the atmosphere is longer by the **airmass** factor,
computed here via the Kasten & Young (1989) empirical formula
(`kasten_young_airmass`), which — unlike the simple `1/cos(zenith)`
plane-parallel approximation — remains well-behaved near the horizon.

---

## 6. Slant-path turbulence & scintillation (`src/atmosphere/turbulence.py`)

**The Hufnagel-Valley profile.** Ground-level Rytov variance
(`rytov_variance` in `src/fso/channel_models.py`) assumes constant `Cn2`
along a horizontal path; a slant path to a satellite instead passes through
turbulence that varies strongly with altitude. The Hufnagel-Valley (HV)
model (`cn2_hv`) is the standard parametric profile for this variation,
combining a high-altitude term driven by wind shear, a mid-altitude term,
and a near-ground term whose strength is set by the single parameter `A0`
(site-dependent ground-level turbulence).

**The Fried parameter.** `fried_r0` computes the Fried parameter `r0` — the
effective aperture diameter over which the atmosphere introduces about one
radian of RMS wavefront phase distortion — via the standard Fried (1966)
path-integral formula through the HV profile, scaled to the slant path by
the secant of the zenith angle. Smaller `r0` means stronger turbulence.

**Slant-path Rytov variance.** The horizontal-path Rytov variance formula
does not directly generalize to a slant path with an altitude-varying
`Cn2(h)`; the correct generalization is a path integral,

```
sigma_R^2 = 2.25 * k^(7/6) * sec(zeta)^(11/6) * int_0^L Cn2(h) * (h-h0)^(5/6) dh
```

for a **plane wave** (`sigma_R2_slant_planewave`). For a ground-to-satellite
**uplink**, however, the correct wave model is **spherical**, not plane —
the beam originates from a point source (the ground telescope) rather than
arriving as a plane wave — and the spherical-wave horizontal-path
coefficient (0.5) differs from the plane-wave coefficient (1.23) by a factor
of `0.5/1.23 ≈ 0.4065`. `sigma_R2_slant` applies this correction
(`SPHERICAL_UPLINK_FACTOR`) to the plane-wave path integral, and is the
function used everywhere downstream in this package, since SNS-TF-QKD here
transmits ground -> satellite. Both the plane-wave path-integral machinery
and the corrected spherical-wave output are validated against the known
constant-`Cn2` horizontal-path closed forms (see `docs/module_documentation.md`
and the smoke tests).

**Provisional PAT quantities.** `aoa_jitter_urad` (Tyler 1994 tracking-error
formula) and `strehl_uncompensated` (Marechal/Noll approximation) give
order-of-magnitude angle-of-arrival jitter and uncompensated Strehl-ratio
estimates from the same HV integral, feeding the receiver-aperture design
discussion in the PAT module.

---

## 7. Seasonal extension (`src/atmosphere/seasonal.py`)

A purely data-driven extension: winter/summer/monsoon aerosol
(AOD/Angstrom-exponent) and ground-turbulence (`A0`) parameter sets per
city, each independently sourced and confidence-tagged, feeding the same
`tau_aerosol_zenith`/`cn2_hv` machinery above with season-specific inputs
rather than the annual-baseline values.

---

## 8. PAT, optics and detector real-loss module (`src/pat/losses.py`)

Replaces the idealized detector-efficiency/dark-count placeholders used in
the fiber/FSO core with a real, itemized hardware loss budget:

**Diffraction / collection efficiency.** The Friis equation for a
point-to-point optical link relates received power to transmitted power
through the product of transmit and receive antenna gains and the
free-space path loss:

```
eta_diffraction = G_tx * G_rx * (lambda / (4*pi*L))^2
```

where `G_tx = 8/theta_div^2` (transmit gain from the beam's far-field
divergence half-angle) and `G_rx = (pi^2/lambda^2)*(D_rx^2 - D_obscured^2)`
(receive gain from the primary aperture area, minus any central
obscuration from a secondary mirror).

**Pointing loss.** Residual pointing jitter (from an imperfect
pointing-acquisition-tracking system) causes the receiver to sit off the
beam's peak intensity some of the time. For a Gaussian beam profile and
Gaussian-distributed jitter, the time-averaged pointing efficiency has a
closed form (`eta_pointing`) in terms of the beam divergence, the RMS
jitter, and any static pointing bias.

**Detector realism.** Real single-photon avalanche photodiodes have a
finite quantum efficiency (`ETA_DET_REAL`), a nonzero intrinsic dark-count
rate, and finite timing jitter/afterpulsing — all itemized from a
manufacturer/publication specification rather than left at the fiber
model's convenient defaults.

**Background light.** A satellite receiver's field of view collects ambient
sky brightness in addition to the signal; this is modeled as an additional
per-pulse "dark count"-like click probability (`R_noise/R_tx`), scaled by a
day/night and per-city multiplier reflecting real urban sky-brightness
differences. Because this term enters the key-rate formula identically to
`p_dark`, it is the dominant real-world determinant of whether a given pass
supports a positive key rate at all (see the root-cause diagnostic in
`src/experiments/pat_experiments.py`).

---

## 9. Coupling to the key-rate engine (`src/coupling/phase6.py`)

Folding turbulence fading (Section 3/6 above) into the fully-dressed
real-loss channel (Section 8) requires ensemble-averaging the key rate over
the Gamma-Gamma fading distribution **with the real, suppressed `p_dark`**
substituted for the fiber model's idealized default — `optimize_ensemble_real`
is this combination: Monte Carlo sampling of the fading channel, joint
`(eps, mu')` optimization at each evaluation, with `p_dark` threaded through
explicitly rather than relying on a module-level default.

**Background-suppression breakeven.** Because background light so strongly
determines link viability, a natural design question is: by what factor
would the background-induced click probability need to be suppressed
(narrower spectral filter, smaller field of view, better spatial/temporal
filtering) for the link to break even (`R=0`)? This is found by a root
search (`scipy.optimize.brentq`) over the suppression factor, at fixed
channel geometry — the standard question a real link-budget design would
ask before specifying hardware.

---

## 10. Weather-gated availability (`src/weather/availability.py`)

A satellite optical link additionally requires **cloud-free line of sight**
— unlike RF links, clouds essentially block the link entirely. This module
converts real station rainy-day-count climatology (day counts, not
intensity, since day-level cloud presence is the relevant proxy for optical
outage) into a monthly availability fraction, with an additional explicit
winter-fog outage term for Delhi (fog is a distinct, non-precipitating
outage mechanism common on the Indo-Gangetic Plain in winter). Combined with
the orbital pass frequency (from Part III's propagator, swept across a
full year) and the independence assumption `joint_availability =
avail_A * avail_B`, this gives the expected number of fully-clear passes
per year and per season.

---

## 11. Joint signal/decoy optimization (`src/optimization/joint.py`)

Two simplifications carried through the fiber/FSO/satellite core so far are
corrected here:

1. **The epsilon bound.** The protocol's send/don't-send sending probability
   `eps` is symmetric under `eps <-> 1-eps` by construction (swapping which
   outcome is labeled "send"), so the physically meaningful domain is `eps
   in (0, 0.5)` — bounding the search at `(0.01, 0.99)`, as the earlier
   optimizer wrappers do, wastes half the search space on a redundant
   mirror image.
2. **The decoy intensity.** Earlier modules fix `mu1 = 0.2*mu2` (a
   conventional, but not optimized, decoy ratio). `sns_key_rate_eta_joint`
   exposes `mu1` as an independent free parameter, and
   `optimize_joint_signal_decoy`/`optimize_joint_ensemble` jointly search
   over `(eps, mu2, mu1)` — three free parameters instead of two — since the
   optimal decoy intensity depends on the specific channel loss and
   background level, not a fixed ratio.

`refactor_correctness_check` validates that the generalized formula, at the
old fixed ratio, reproduces the original formula's output exactly — the
correct check for a genuine generalization rather than a subtly different
new formula.

---

## 12. Multi-city links and the asymmetric protocol (`src/links/`)

**Why per-station rates and a simple product are both wrong.** A first
instinct for evaluating a link between two *different* cities (with
different channel transmittances `eta_A != eta_B`, since SNS-TF-QKD here is
a joint three-party protocol) might be to (a) treat each ground station as
having its own independent rate, or (b) multiply the two arms' transmittances
into one effective `eta = eta_A * eta_B` and reuse the symmetric formula.
Both are physically wrong:

- **(a)** is wrong because Alice and Bob do not have separate rates — they
  jointly feed Charlie's *one* shared interference measurement; there is no
  such thing as "Alice's key rate" independent of Bob's channel.
- **(b)** is wrong for a deeper reason: twin-field protocols get their
  square-root-of-total-loss scaling advantage over ordinary repeaterless QKD
  precisely because each arm contributes **linearly** (through the combined
  optical field amplitude) to Charlie's click probability, not
  multiplicatively. Multiplying the two arms' transmittances together first
  and treating the product as a single effective `eta` would silently erase
  that linear structure and collapse the protocol back to ordinary,
  worse-scaling linear-loss behavior.

**The correct generalization** (`src/links/metrics.py`) keeps each arm's
transmittance separate through every stage of the rate formulas:

- The **X-window** gain combines both arms additively in the click-rate
  exponent (`mu*eta_A + mu*eta_B`), since both parties always send in the
  X-window and either arm's photon can trigger Charlie's detector.
- The **Z-window** gain and QBER are split into the four true joint cases —
  Alice-only-sends, Bob-only-sends, both-send, neither-sends — each with its
  own click probability, rather than averaged or multiplied together, since
  only one of Alice/Bob independently deciding to send is the physically
  correct picture of the protocol's own sending-probability mechanism.
- Background light is **summed** across both ground stations' own collected
  stray light before reaching Charlie's one shared detector, while Charlie's
  dark-count rate stays single (one physical measurement system).

`validate_symmetric_limit` confirms this generalization reduces exactly to
the original symmetric formulas when `eta_A = eta_B`, which is the correct
correctness criterion — a genuine generalization must contain the original
special case exactly, not merely produce plausible-looking new numbers.

**Independent per-arm fading.** `optimize_joint_asym_ensemble` draws
*independent* Gamma-Gamma fading samples for each arm (Delhi and Mumbai, for
example, sit under different turbulence conditions and should not share one
averaged fading process), then jointly optimizes `(eps, mu2, mu1)` under
this two-arm ensemble average — the natural generalization of Section 11's
single-arm ensemble optimizer.

**Per-link background suppression.** Because different links have different
combined background levels and different channel losses, the
background-suppression breakeven factor (Section 9) is re-derived
independently for each of the four evaluated city-pair links rather than
reused from the Delhi-Mumbai case, with a uniform design margin (currently
3x) adopted above breakeven.
