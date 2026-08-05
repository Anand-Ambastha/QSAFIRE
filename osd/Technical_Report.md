# Sending-or-Not-Sending Twin-Field Quantum Key Distribution: Complete System Design and Simulation Architecture Report

**A Technical Documentation for Full System Reconstruction in OptiSystem/MATLAB/Python**

*Prepared as laboratory-grade technical documentation*

---

## Preface and Scope of This Document

This report documents, component by component, the OptiSystem simulation architecture implementing a Sending-or-Not-Sending Twin-Field Quantum Key Distribution (SNS-TFQKD) link between two transmitting stations, Alice and Bob, and a central untrusted measurement station, Charlie. The document is written for a reader who already possesses a working knowledge of optical communication theory and quantum-optical concepts, but who has never seen this particular architecture and must be able to rebuild it from scratch — including every parameter, every block, every MATLAB co-simulation stage, every CSV artifact, and every hand-off to the external Python post-processing environment.

The document intentionally separates concerns along the same boundary that the simulation itself enforces: OptiSystem and its embedded MATLAB co-simulation blocks are responsible exclusively for physical-layer signal generation, propagation, and detection; nothing resembling key distillation, error-rate estimation, or secure-key-rate computation is performed inside the simulation. All protocol-layer analytics are deferred to an external Python environment that ingests the CSV artifacts produced during the OptiSystem run. This separation is not incidental — it is the central architectural decision of the whole design, and it is revisited repeatedly throughout the chapters below because it explains why certain blocks exist, why certain data are logged rather than acted upon, and why no decision block in the schematic ever "closes the loop" on key generation.

Wherever MATLAB source code would ordinarily appear, this report inserts a labelled placeholder — *MATLAB Implementation* — instead of code, and describes in prose, pseudocode, and flowchart form what that code must do. This mirrors the constraint under which this documentation was produced and keeps the report usable as a specification document independent of a particular MATLAB coding style.

---

# Chapter 1 — Introduction

## 1.1 Background of Quantum Key Distribution

Quantum key distribution (QKD) allows two parties, conventionally named Alice and Bob, to establish a shared secret key whose secrecy is guaranteed not by computational hardness assumptions but by the physical laws governing quantum systems. The founding protocol, BB84, encodes classical bits onto the polarization or phase of single photons prepared in one of two mutually unbiased bases; any attempt by an eavesdropper (Eve) to gain information about the transmitted bits necessarily perturbs the quantum state in a way that is, in principle, statistically detectable by Alice and Bob during a public reconciliation phase. This "observation disturbs the system" property, rooted in the no-cloning theorem and in the Heisenberg uncertainty relation, is what elevates QKD from a merely clever cryptographic protocol to one with an information-theoretic security proof, provided the physical implementation matches the theoretical assumptions closely enough.

In practice, QKD is implemented over optical fibre or free-space optical links, using weak coherent pulses (WCPs) emitted by attenuated lasers as a practical stand-in for genuine single-photon sources. This substitution reintroduces a vulnerability — the photon-number-splitting (PNS) attack — because a coherent state occasionally contains two or more photons, and Eve can, in principle, siphon off the surplus photons without perturbing the single-photon component that Alice and Bob believe they are exchanging. The decoy-state method neutralises the PNS attack by having Alice randomly vary the mean photon number of her pulses among several intensity levels (signal, decoy, vacuum) so that Bob's detection statistics as a function of intensity betray any photon-number-dependent attack strategy, allowing Alice and Bob to place a tight statistical bound on the single-photon-state fraction and the associated error rate.

Even with decoy states solving the PNS problem, standard point-to-point QKD is subject to a hard, information-theoretic ceiling on achievable secret-key rate as a function of channel loss: the secret-key capacity (SKC) of a lossy bosonic channel, which scales *linearly* with the channel transmittance η. Because η decays exponentially with fibre length (η = 10^(−αL/10) for attenuation coefficient α, typically 0.2 dB/km for standard single-mode fibre), the achievable key rate of any repeaterless point-to-point QKD scheme collapses exponentially with distance. This is the rate–distance limit, and prior to 2018 it was believed that surpassing it required a genuine quantum repeater — a device capable of quantum memory, entanglement swapping, and/or quantum error correction, none of which is presently deployable outside laboratory demonstrations.

## 1.2 Evolution from BB84 to TF-QKD

The architecture documented in this report descends directly from the twin-field QKD (TF-QKD) proposal of Lucamarini, Yuan, Dynes and Shields (*Nature* 557, 400–403, 2018). The core insight of TF-QKD is deceptively simple: instead of Alice sending a photon all the way to Bob (distance L), or Alice and Bob each sending a photon to a possibly-untrusted intermediary Charlie who performs a Bell-state-like coincidence measurement (as in measurement-device-independent QKD, MDI-QKD), TF-QKD has Alice and Bob each transmit *phase-coherent, phase-randomized weak coherent pulses* to Charlie, positioned midway between them, where the two fields are made to interfere on a single beam splitter. A single-photon detection event on one of Charlie's two output ports (rather than a two-photon coincidence, as required by MDI-QKD) heralds an effective quantum-correlated event between Alice's and Bob's encoded bits.

Because the detection event now depends on the interference of two *fields*, each of which need only survive a one-way trip of length L to Charlie (rather than 2L point-to-point, or the two-photon joint probability of MDI-QKD), the click rate — and therefore the key rate — scales with the *square root* of the total channel transmittance rather than linearly with it. This single change in scaling law is the reason TF-QKD is able to beat the repeaterless SKC bound using only equipment that already exists in decoy-state MDI-QKD implementations: lasers, intensity modulators, phase modulators, variable optical attenuators, single-mode fibre, a 50:50 beam splitter, and single-photon detectors.

The price paid for this square-root scaling is that Alice's and Bob's optical fields must remain mutually phase-coherent (to within a controllable global-phase slice) over the entire propagation distance, which in the original proposal is enforced by transmitting an interleaved sequence of bright reference pulses (used for active phase-locking and drift compensation at Charlie's station) and dim quantum-regime pulses (used for the actual key-encoding interference). The original paper further shows experimentally that phase drift over hundreds of kilometres of standard fibre, though non-negligible, is slow enough (a few rad/ms) to be tracked and compensated with a feedback loop acting on Charlie's local phase modulator, and that interference visibility in excess of 99.6% is achievable even at 550 km.

## 1.3 Need for the SNS Protocol

The original TF-QKD paper explicitly flagged an unresolved security gap: because Alice's and Bob's bit values are encoded on the *global random phase* ρ_a, ρ_b of their coherent pulses, and because that global phase must eventually be revealed (in coarse "phase-slice" form) during the public sifting stage so that Alice and Bob can identify which pulse pairs were mutually phase-compatible ("twins"), an eavesdropper who exploits the *timing* of that phase announcement can, in specific attack constructions, extract information about the encoded bit without introducing a detectable disturbance. Wang, Yu and Hu (arXiv:1805.09222) constructed an explicit non-destructive, sequential eavesdropping attack of exactly this kind and showed that the traditional decoy-state key-rate formula, applied naively to the original TF-QKD protocol, overestimates the achievable secure key rate — in the worst case it returns a nonzero key rate for a channel on which Eve in fact possesses full information.

The Sending-or-Not-Sending (SNS) protocol, introduced in the same paper, closes this gap by a structural change rather than a patch: in the Z-basis (the basis actually used to generate key bits), Alice and Bob no longer encode a bit on the value of a randomly chosen global phase that must later be disclosed. Instead, each of them makes a *local, private, binary decision* — "send" a weak coherent pulse of fixed intensity μ′, or "do not send" (i.e., send vacuum) — with sending probability ε. Because the phase information of Z-basis (signal) pulses is *never* publicly announced, the decoy-state analysis that Wang et al. use to bound the single-photon yield and phase-error rate applies without modification, and the specific eavesdropping construction that broke the original protocol's security proof no longer has a foothold. A separate, disjoint set of time windows — the X-basis or "decoy" windows — is reserved purely for phase-randomized decoy pulses used to estimate the single-photon interference misalignment error, and it is only in these X-windows that the delicate long-distance single-photon interference must succeed with reasonable visibility.

This separation of duties between the Z-basis (key generation, no interference requirement, negligible intrinsic error) and the X-basis (parameter estimation, interference required, error-tolerant) is the reason the SNS protocol tolerates single-photon misalignment error rates up to and beyond 45%, whereas the original TF-QKD proposal required interference visibility good enough to keep its intrinsic QBER at only a few percent. This is a decisive practical advantage: phase-tracking hardware that would be marginal for standard TF-QKD becomes comfortably adequate for SNS-TFQKD, which is precisely why SNS-TFQKD, rather than the original TF-QKD protocol, is the variant implemented in the OptiSystem architecture documented in this report.

## 1.4 Measurement-Device-Independent QKD

Both TF-QKD and SNS-TFQKD inherit the measurement-device-independent (MDI) property first introduced by Lo, Curty and Qi in 2012: neither Alice nor Bob need trust Charlie or his detectors. Because the secret bit is inferred solely from *which* of Charlie's two detectors clicked (and, in the SNS case, from Alice's and Bob's own private sending/not-sending records), any detector-side attack — blinding attacks, detector-efficiency-mismatch attacks, time-shift attacks, and so on — that would compromise a conventional QKD receiver is automatically neutralised, because Charlie's detectors are assumed to be entirely under an adversary's control in the security proof. This is architecturally significant for the design documented here: the single-photon detectors (SPDs) at Charlie's station are treated, throughout, purely as *click/no-click event loggers*. No decision logic anywhere downstream of the SPDs trusts the detectors' internal calibration; all detector imperfections (dark counts, finite efficiency) are instead folded into the statistical model that the external Python post-processing stage uses to compute bounds on the secure key rate.

## 1.5 Objectives of This Architecture

The OptiSystem architecture documented here has the following explicit objectives:

- Simulate, at the physical layer, two independent SNS-TFQKD transmitter chains (Alice, Bob) operating at 10 Gbit/s symbol rate, each producing decoy and signal pulses according to the SNS protocol's Z-window/X-window structure.
- Propagate both transmitters' optical output through independent, identical 40 km spans of standard single-mode fibre (SMF-28) toward a central interference station, Charlie.
- Perform the interfering beam-splitter measurement at Charlie's station using two arms (Z-window and X-window processing paths, each duplicated for redundancy/parallel channel handling as shown in the schematic), each terminated with a pair of single-photon detectors.
- Recover the detection events, correlate them against transmitted bit sequences via dedicated MATLAB and Data Recovery blocks, and log every relevant classical variable (bit index, window type, sent/not-sent decision, intensity class, random phase, detector click pattern, arrival time) into synchronized CSV files.
- Explicitly avoid performing any protocol-layer computation (sifting, gain, QBER, secure key rate, finite-key statistics) inside OptiSystem, deferring all such computation to an external Python pipeline that consumes the CSV artifacts.

## 1.6 Applications

SNS-TFQKD architectures of the kind modelled here are directly relevant to: metropolitan and inter-city quantum-secure backbone links where a trusted, centrally located relay node (Charlie) is operationally convenient; quantum-secured data-centre interconnects; QKD network topologies where an untrusted or semi-trusted central switching node is unavoidable for economic reasons (single fibre plant serving multiple access points); and testbeds for evaluating detector technologies and phase-stabilization schemes intended for eventual field deployment over hundreds of kilometres of deployed fibre.

## 1.7 Advantages

The principal engineering advantages of the SNS-TFQKD architecture, relative to both conventional point-to-point decoy-state BB84 and to decoy-state MDI-QKD, are: (i) square-root-of-transmittance key-rate scaling, extending the maximum secure distance considerably beyond the repeaterless SKC bound; (ii) full measurement-device independence, removing an entire class of detector side-channel attacks; (iii), specific to the SNS variant, an intrinsically negligible Z-basis error rate because no single-photon interference is required for key-bit generation, decoupling the achievable secure distance from the technologically demanding requirement of maintaining sub-percent-level interference visibility over hundreds of kilometres; and (iv) full compatibility with the standard, mature decoy-state analysis toolchain, because the phase information of signal (Z-basis) pulses is never disclosed.

## 1.8 Limitations

The architecture also carries limitations that a rebuilder must understand: the protocol still fundamentally relies on **weak coherent pulses**, so photon statistics remain Poissonian and decoy-state estimation, while tight, is not identical to a true single-photon source; the **X-basis** windows still require single-photon-level interference and therefore still demand active phase tracking and compensation infrastructure at Charlie's station, even though the tolerable misalignment error is much larger than in the original TF-QKD scheme; the **security proof** of the SNS protocol (Section V of Wang et al.) is intricate, resting on a virtual-protocol / tagging-model reduction, and any deviation of the real implementation from the assumed statistical model (e.g., correlated phase noise between Z- and X-windows, imperfect intensity modulator extinction ratio, detector afterpulsing) can erode the margin between the simulated and the provably secure key rate; and, at the level of this specific OptiSystem model, the **finite fibre length (40 km per span, four spans total)** used in the simulation is a laboratory-scale placeholder chosen for tractable simulation run-time, not a proxy for the multi-hundred-kilometre distances at which SNS-TFQKD's real advantage over MDI-QKD becomes apparent; a rebuilder targeting realistic long-haul performance must substitute longer SMF-28 spans and correspondingly adjust detector dark-count and dead-time parameters.

---

# Chapter 2 — Complete System Architecture

## 2.1 Overview of the Signal Flow

The system is organized into three physical stations — **Alice**, **Bob**, and **Charlie** — connected by four independent spans of single-mode fibre (SMF-28, 40 km each), and it is built from four cooperating computational domains: the **electrical domain** (bit-sequence generation and the classical decision logic implemented as MATLAB co-simulation components), the **optical domain** (laser sources, modulators, fibre propagation, beam-splitting interference, and single-photon detection), the **MATLAB processing layer** (which sits between the electrical and optical domains at multiple points, injecting protocol logic that OptiSystem's native electrical/optical block library cannot express), and the **offline Python layer** (which performs everything downstream of raw detection-event logging).

At the top level, each of Alice and Bob instantiates two parallel processing branches, labelled **Z-WINDOW** and **X-WINDOW** in the schematic, corresponding exactly to the SNS protocol's Z-basis (signal, key-generating) and X-basis (decoy, parameter-estimating) time windows. Both branches share a single upstream bit-sequence generator and a single Window Selection MATLAB component that partitions the simulation's bit stream into Z-window and X-window sub-populations. Each branch then independently modulates a continuous-wave (CW) laser source: the Z-window branch applies SNS Logic (the send/not-send decision), while the X-window branch applies Decoy Estimation (choice among signal/weak-decoy/vacuum intensities) and Random Phase modulation (uniform phase randomization required for X-basis interference-based parameter estimation). Both branches recombine optically before entering the shared 40 km fibre span connecting each user to Charlie.

At Charlie's station, the two incoming optical fields (one from Alice, one from Bob) are combined at an interfering (X) coupler acting as the 50:50 beam-splitter of the TF-QKD/SNS scheme, whose two output ports are each monitored by a pair of single-photon detectors (SPD) — one pair dedicated to the Z-window processing path, the other pair to the X-window processing path, mirroring the transmitter-side branch structure. Each detector's click stream is fed into a Data Recovery block (which reconstructs the timing and logical click/no-click sequence), visualized on Dual Port Oscilloscope and Dual Port Binary Sequence Visualizers for engineering-verification purposes, and ultimately logged, together with all locally available window/intensity/phase metadata, into the final synchronized CSV artifact via a dedicated MATLAB logging component (CSV FINAL Z-WINDOW – MATLAB, and its analogous instance for the X-window path).

## 2.2 Electrical Domain

The electrical domain is responsible for everything that happens before light leaves a laser: pseudo-random and user-defined bit-sequence generation, NRZ pulse shaping of those bit sequences into electrical drive waveforms, and the classical decision logic (window selection, SNS sending decision, decoy-intensity class selection, random-phase index selection) that OptiSystem's MATLAB co-simulation component executes numerically, sample-by-sample, in lock-step with the optical simulation clock. Every electrical waveform in this architecture is a rectangular NRZ pulse train at the system bit rate (10 Gbit/s); no pulse shaping filter (raised-cosine, Bessel, Gaussian) is applied, because the electrical waveforms in this design are used purely as *digital control signals* driving optical modulators in an on/off or discretely-stepped fashion, not as analog information-bearing waveforms whose spectral shape must be preserved.

## 2.3 Optical Domain

The optical domain comprises: five independent CW laser sources (one for Alice's Z-window branch, one for Alice's X-window branch, one for Bob's Z-window branch, one for Bob's X-window branch — noting that the schematic instantiates "CW Laser", "CW Laser_1", "CW Laser_15", and "CW Laser_16" — plus any additional laser instances required for full branch independence), each operating at the standard telecom wavelength corresponding to 193.1 THz (≈1552.52 nm, ITU-T C-band channel) and at 0 dBm optical output power; a Mach–Zehnder Modulator (Analytical model) on each Z-window branch, used to imprint the NRZ-shaped bit pattern as an on/off amplitude pattern; Intensity Modulators on both Z- and X-window branches, used respectively to implement the SNS Logic's binary send/not-send gating (Z-window) and the Decoy Estimation module's multi-level intensity selection (X-window); Phase Modulators on the X-window branches, driven by the Random Phase MATLAB component, to imprint the uniformly distributed global phase required for X-basis interference statistics; four spans of SMF-28 fibre (40 km each) carrying, respectively, Alice's combined optical output to Charlie, and Bob's combined optical output to Charlie, with each user's Z-window and X-window optical paths kept on physically distinct fibre runs in this model (as evidenced by the four separate SMF-28 blocks in the schematic) rather than being wavelength- or time-multiplexed onto a single fibre; two X Couplers at Charlie's station, each acting as the interfering beam splitter for one basis (Z or X); and four Single Photon Detectors (SPD, SPD_1, SPD_2, SPD_3), two per basis, monitoring the two output ports of each X Coupler.

## 2.4 MATLAB Processing

MATLAB co-simulation components are the mechanism by which this architecture expresses protocol logic that has no native representation in OptiSystem's built-in block library. OptiSystem is fundamentally a physical-layer, waveform-level simulator; it has no native concept of "with probability ε, do not transmit this symbol" conditioned on an upstream random bit, nor any native concept of "select one of three intensity classes according to a pre-agreed decoy-state probability distribution and log the choice to a CSV row indexed consistently with five other MATLAB blocks running elsewhere in the schematic." Each MATLAB Component in this design therefore acts as a *programmable, stateful, arbitrary-logic tap* inserted into the OptiSystem signal-flow graph: it accepts one or more OptiSystem signals (electrical and/or optical) as MATLAB input variables, executes ordinary MATLAB code against them (including file I/O for CSV logging, which is why these components are able to accumulate a persistent log across the whole simulation run rather than being limited to per-sample processing), and returns one or more signals of the appropriate OptiSystem signal type back into the simulation graph.

Six named MATLAB component roles appear in this architecture, several instantiated twice (once for Alice, once for Bob): **Window Selection MATLAB** (routes each bit index into the Z-window or X-window branch based on the PRBS bit), **SNS Logic MATLAB** (implements the send/not-send decision within Z-windows), **Decoy Estimation MATLAB** (implements the signal/weak-decoy/vacuum intensity selection within X-windows), **Random Phase MATLAB** (generates the uniformly distributed global phase for X-window pulses), and **CSV FINAL [Z/X]-WINDOW MATLAB** (the terminal logging stage at Charlie's station that consolidates detector click data with the upstream transmitted-side metadata into the final synchronized CSV files). Every one of these is documented in full in Chapters 4 through 6 below, including inputs, outputs, internal algorithm, pseudocode, and CSV schema; consistent with the constraints of this report, no MATLAB source code is reproduced — only placeholders marking where it must be inserted.

## 2.5 Detection

Detection in this architecture is deliberately "dumb" at the OptiSystem level: the SPD components report only click/no-click events (with associated arrival-time and, where modelled, dark-count statistics) on each of Charlie's four detector ports. No logic anywhere in the OptiSystem schematic decides, on the fly, whether a given detection event constitutes a valid SNS "effective event," nor does any block compute a bit value from it. That interpretation step — matching a detector click against the SNS protocol's effective-event criteria (Section II, Step 2 of Wang et al.) — is explicitly deferred to the offline Python stage, which has simultaneous access to Alice's CSV, Bob's CSV, and Charlie's detector CSV and can therefore perform the necessary three-way join.

## 2.6 Synchronization

Every logged event in this architecture, whether it originates in Alice's SNS Logic MATLAB block, Bob's Decoy Estimation MATLAB block, or Charlie's detector logging stage, is indexed by a common, monotonically increasing **bit index** (equivalently, time-window index) that is derived from the same underlying simulation clock (bit rate = 1×10^10 bit/s, samples per bit = 32, giving a fixed sample-domain relationship between simulation time and bit index for every block in the schematic). Because all MATLAB components execute in lock-step with the OptiSystem sample clock and because none of them independently re-samples or re-clocks the bit stream, the bit-index column functions as a global synchronization key across every CSV file the architecture produces; Chapter 7 develops this in full detail, including the failure modes (dropped samples, off-by-one indexing errors introduced by guard bits) that a rebuilder must guard against.

## 2.7 Data Logging and CSV Generation

Each station accumulates its own consolidated CSV file over the course of the simulation run: **Alice.csv**, **Bob.csv**, and one or more detector-side CSV files at Charlie's station (conceptually **Detector.csv**, materialized in the schematic as the CSV FINAL Z-WINDOW MATLAB output and its X-window analogue). These files are written incrementally, one row appended per processed time window, by MATLAB's ordinary file-I/O functions operating in "append" mode from inside each MATLAB Component's per-sample or per-symbol callback. Because OptiSystem re-invokes a MATLAB Component's code once per relevant signal sample (or once per symbol period, depending on how the component's timing is configured), the CSV-writing code must be written to append exactly one row per meaningful bit index, not one row per raw waveform sample — this is one of the most consequential implementation details in the entire architecture and is discussed at length in Chapter 4.

## 2.8 Offline Processing

Offline processing is, by design, everything that happens after the OptiSystem run completes and the CSV files are closed. This includes CSV loading and validation, cross-file synchronization and index alignment, sifting (discarding time windows where Alice and Bob's basis choices did not match, or where Charlie did not report a valid single-click effective event), decoy-state yield and error-rate estimation (implementing, e.g., Equations 44–45 of Wang et al.), gain and QBER computation, error-correction-cost accounting, finite-key statistical bounding, and the final asymptotic or finite-key secure-key-rate figure (Equations 3–4 of Wang et al.). None of this logic exists inside the OptiSystem schematic; Chapter 8 documents the full offline pipeline in detail, framed explicitly as a specification for the external Python implementation rather than as a description of any block visible in the OptiSystem canvas.

## 2.9 Why OptiSystem and MATLAB Are Integrated

OptiSystem is chosen as the physical-layer engine because it natively models the propagation-level physics that matter for this system — laser linewidth and RIN, modulator extinction ratio and chirp, fibre attenuation and dispersion, detector quantum efficiency, dark-count rate and jitter — using validated, parameterized component models rather than requiring the rebuilder to hand-code optical propagation from first principles. MATLAB co-simulation is integrated at every point where the SNS protocol's logic exceeds what a linear, waveform-oriented simulator can express natively: probabilistic branching per time window, persistent state across the whole simulation run (e.g., an evolving CSV log, or a decoy-intensity sequence drawn from a pre-agreed probability table), and the injection of externally generated random numbers (phase values, sending decisions) that must be logged for later, deterministic replay/audit by the offline analysis stage. The result is a division of labour in which OptiSystem "owns" continuous-time optical and electrical waveform physics, while MATLAB "owns" discrete-time, per-window protocol decisions and their bookkeeping.

## 2.10 Why Python Is Used Only for Final Post-Processing

Python is deliberately excluded from any role inside the real-time OptiSystem/MATLAB co-simulation loop. There are three reasons for this separation. First, performance: OptiSystem's MATLAB co-simulation interface is optimized for tight, low-latency round trips with MATLAB specifically; introducing a second, cross-process hand-off to a Python interpreter on every simulated sample would be prohibitively slow for a run producing tens of thousands of samples. Second, architectural clarity: keeping all protocol-layer cryptographic analysis (gain, QBER, secure-key-rate) *outside* the physical-layer simulator makes it straightforward to verify that the physical-layer model is not, even inadvertently, "peeking" at information that a real Charlie, Alice, or Bob would not have access to during the actual quantum transmission phase — a property that matters for the model's fidelity as a stand-in for a real experimental or field system. Third, reusability: because the CSV artifacts fully capture every random draw and every detection event, the same CSV output can be replayed through multiple independent Python analysis pipelines (e.g., an asymptotic key-rate calculator and a separate finite-key statistical analysis tool) without ever re-running the OptiSystem simulation, which is by far the most computationally expensive stage of the whole workflow.

---

# Chapter 3 — Simulation Parameters

## 3.1 Global Simulation Layout Parameters

The table below documents every global parameter visible in the OptiSystem layout header, together with its meaning, units, the rationale for the chosen value, its impact on the simulation, and the principal trade-off a rebuilder must weigh when changing it.

| Parameter | Value | Units | Meaning |
|---|---|---|---|
| Bit rate | 1 × 10¹⁰ | bit/s | Symbol clock driving every PRBS, NRZ, and modulator block in the schematic |
| Sequence length | 1024 | bits | Total number of bits generated per PRBS/user-defined sequence per simulation run |
| Samples per bit | 32 | samples/bit | Oversampling factor used to represent each bit-period waveform |
| Sample rate | 3.2 × 10¹¹ | Hz | Product of bit rate and samples per bit; the master waveform sampling clock |
| Number of samples | 32768 | samples | Sequence length × samples per bit; total waveform vector length per run |
| Symbol rate | 1 × 10¹⁰ | symbols/s | Equal to bit rate, since each transmitted symbol here is a single binary decision |
| Time window | 1.024 × 10⁻⁷ | s | Total simulated time span = sequence length / bit rate |
| Guard bits | 0 | bits | Number of bits excluded from analysis at sequence boundaries to avoid filter transients |

### 3.1.1 Bit Rate — Meaning, Reason for Selection, Impact, Trade-offs

The bit rate of 10 Gbit/s sets the fundamental time-window duration for every SNS protocol decision in the model: each 100 ps interval corresponds to exactly one "time window i" in the language of Wang et al.'s protocol description, i.e., one opportunity for Alice and Bob to each independently choose a Z-window (signal) or X-window (decoy) commitment and, within a Z-window, to decide send or not-send. This value is chosen because it corresponds to a realistic, currently deployable clock rate for high-speed intensity- and phase-modulator drive electronics and single-photon-detector gating circuitry, making the simulation representative of a physically buildable system rather than an idealized abstraction. Increasing the bit rate directly increases the raw pulse-emission rate and therefore the achievable raw key rate (all else equal), but it also proportionally tightens the timing budget available for phase-modulator settling, detector gating, and — critically — Charlie's active phase-drift compensation loop, which must track and correct phase noise fast enough relative to the pulse repetition period. Decreasing the bit rate relaxes these timing constraints at a direct, linear cost in raw key rate.

### 3.1.2 Sequence Length and Number of Samples

A sequence length of 1024 bits, at 32 samples per bit, yields exactly 32768 waveform samples — a power-of-two length that is favourable both for FFT-based operations OptiSystem may perform internally (chromatic dispersion modelling in the fibre block, for instance, typically uses frequency-domain propagation) and for straightforward indexing inside the MATLAB co-simulation components. This sequence length is a *simulation-scale* choice, not a protocol-scale one: 1024 raw bits are far too few to obtain statistically meaningful decoy-state bounds or finite-key security margins, so a rebuilder intending to reproduce publication-quality asymptotic key-rate curves (of the kind shown in Wang et al.'s Figures 1–2) must either increase this value by orders of magnitude or, more practically, run the OptiSystem simulation repeatedly (via the "Sweep Iteration" mechanism referenced in the layout header) and concatenate the resulting CSV logs before handing them to the offline Python statistical pipeline. Increasing the sequence length increases simulation run-time roughly linearly (and, for MATLAB components performing per-sample CSV writes, can increase run-time super-linearly if the CSV file is re-opened rather than kept open across the run — an implementation detail flagged again in Chapter 4).

### 3.1.3 Samples per Bit and Sample Rate

Thirty-two samples per bit, at the 10 Gbit/s bit rate, fixes the sample rate at 3.2 × 10¹¹ Hz (320 GSa/s). This oversampling factor is chosen to give OptiSystem's optical and electrical filters (and, in particular, the Mach–Zehnder modulator's frequency-domain analytical model) enough resolution to represent the modulated optical spectrum without aliasing, while remaining computationally tractable — doubling the oversampling factor quadruples memory footprint for a fixed sequence length (more samples per bit *and* correspondingly more total samples) and increases run-time correspondingly. Thirty-two samples per bit is a conservative-but-standard choice for NRZ-modulated on/off-type signals of the kind used throughout this architecture (as opposed to more demanding formats such as high-order QAM, which typically require finer oversampling for accurate simulation).

### 3.1.4 Time Window

The 1.024 × 10⁻⁷ s (102.4 ns) total time window is simply sequence length divided by bit rate; it represents the total simulated duration of one OptiSystem run and therefore the maximum time span over which any single MATLAB component's persistent internal state (e.g., a running tally of sent/not-sent decisions, or an open CSV file handle) needs to remain valid within one execution of the simulation. Given the 40 km fibre spans in this design (round-trip-equivalent one-way propagation delay of roughly 196 μs at ~204,000 km/s group velocity in SMF-28), a rebuilder should note that the simulated time window is far shorter than the one-way propagation delay across a single fibre span — meaning that, in a single execution of this OptiSystem layout, symbols transmitted by Alice within the simulated window have not yet, within simulated time, actually arrived at Charlie by the time the run completes its bookkeeping in the naive sense. In practice, OptiSystem's fibre propagation block accounts for this delay internally (it is a linear channel model applying the appropriate group delay to the waveform), so downstream blocks receive correctly time-shifted waveforms; but a rebuilder implementing custom MATLAB timing logic that assumes "same simulated time index equals same physical event" across Alice, Bob, and Charlie must explicitly compensate for the fibre's group delay when correlating each station's log by real time rather than relying purely on the bit-index synchronization scheme described in Chapter 7.

### 3.1.5 Guard Bits

Guard bits are set to zero in this layout, meaning no bits at the start or end of the sequence are excluded from protocol logic or CSV logging. This is an appropriate choice for a system, like this one, in which every block operating on the electrical/logical bit stream is a MATLAB co-simulation component operating symbol-synchronously rather than a linear filter with a transient response (such as a pulse-shaping FIR filter, which would otherwise corrupt the first and last few symbols and necessitate excluding them via a nonzero guard-bit count).

### 3.1.6 Laser Frequency and Power

Every CW laser source in this architecture (CW Laser, CW Laser_1, CW Laser_15, CW Laser_16) is set to 193.1 THz — the ITU-T DWDM grid's canonical 100 GHz-spaced C-band reference channel, corresponding to a vacuum wavelength of approximately 1552.52 nm — and to 0 dBm (1 mW) output power. The 193.1 THz choice reflects standard telecom practice and, crucially, ensures the simulated wavelength falls within the low-attenuation window of SMF-28 fibre (attenuation minimum near 1550 nm, typically ≈0.2 dB/km, consistent with the attenuation coefficient used in the Lucamarini et al. bound calculations). The 0 dBm laser output is a *pre-attenuation* reference power: this is the power at the laser's output facet, well above the single-photon regime; it is the downstream Variable Optical Attenuators and Intensity Modulators (configured as part of the SNS Logic and Decoy Estimation stages) that reduce the mean photon number per pulse down to the μ′ (signal), μ (weak decoy), and vacuum levels required by the protocol. Keeping the laser itself at a conventional, high, easily stabilized power and performing all intensity discrimination downstream is both physically realistic (real diode lasers are not directly modulated down to single-photon output powers with adequate extinction ratio) and numerically convenient in OptiSystem, which handles power quantities more robustly in the classical-power regime before the final attenuation stage.

### 3.1.7 Fibre Length

All four SMF-28 fibre blocks in this architecture are configured to 40 km. This value is a laboratory-scale placeholder selected for simulation tractability (shorter fibre spans keep chromatic-dispersion and phase-noise accumulation small enough that the interference visibility at Charlie's station remains high without requiring an elaborate active phase-compensation MATLAB model) rather than a value intended to demonstrate SNS-TFQKD's headline long-distance performance. A rebuilder targeting the multi-hundred-kilometre regime discussed in both source papers must (a) increase each SMF-28 block's length parameter, (b) verify that OptiSystem's fibre attenuation coefficient is explicitly set to 0.2 dB/km (matching the assumption used throughout the SNS numerical simulations), and (c) introduce or strengthen an active phase-drift compensation stage at Charlie's Phase Modulator, since the phase-drift-rate figures reported experimentally by Lucamarini et al. (2.4 rad/ms at 100 km, 6.0 rad/ms at 550 km) become increasingly significant relative to the 100 ps symbol period as distance grows.

### 3.1.8 Detector Efficiency and Attenuator Values

Although not shown explicitly in the layout header table, the SPD blocks' quantum efficiency and dark-count-probability parameters, and the Variable Optical Attenuator / Intensity Modulator extinction-ratio settings that realize the μ′/μ/vacuum intensity classes, are the parameters that most directly determine the simulated secure-key rate once the CSV logs reach the Python post-processing stage. Consistent with the numerical-simulation parameters used in both source papers (η_det ≈ 30% in the original TF-QKD paper's realistic scenario; 80% detection efficiency and 10⁻¹¹ dark-count probability per pulse in the SNS paper's own numerical simulation), a rebuilder should treat these as the primary "dial" for exploring the trade-off between detector cost/complexity (higher-efficiency superconducting nanowire single-photon detectors versus lower-cost, lower-efficiency avalanche photodiodes) and achievable secure distance, and should log the exact values used in each simulation sweep iteration into the CSV metadata so the offline Python analysis can correctly parameterize its yield and QBER models.

## 3.2 Tradeoff Summary Table

| Parameter increased | Primary benefit | Primary cost |
|---|---|---|
| Bit rate | Higher raw key rate | Tighter phase-tracking and detector-gating timing budget |
| Sequence length | Better decoy-state statistics, lower finite-key penalty | Longer run-time, larger CSV files |
| Samples per bit | Higher waveform fidelity, less aliasing | Higher memory and CPU cost |
| Fibre length | More representative of long-haul deployment | Lower interference visibility, higher phase-drift error, lower raw click rate |
| Laser power (pre-attenuation) | Better SNR for classical/bright-pulse phase-locking channel | No direct effect on quantum-regime pulses once attenuated to μ, μ′ |
| Detector efficiency | Higher click rate, longer achievable secure distance | Increased detector cost/complexity, potentially higher dark counts (technology-dependent) |

---

# Chapter 4 — Alice Transmitter

Alice's transmitter chain begins with bit-sequence generation and terminates at the launch facet of the SMF-28 fibre span connecting her station to Charlie. Every component is documented below with its purpose, inputs, outputs, internal working, relevant theory, governing equations, key parameters and their effects, advantages, limitations, and its interaction with adjacent blocks.

## 4.1 Pseudo-Random Bit Sequence Generator

**Purpose.** Produces the master pseudo-random bit stream from which Alice's Window Selection decision (Z-window versus X-window) is derived for every time window of the simulation.

**Inputs.** None (internally seeded); configured with Bit rate = Bit rate (bit/s), i.e., bound to the global 10 Gbit/s system clock.

**Outputs.** A single electrical bit-sequence signal, one bit per time window, delivered synchronously to the downstream NRZ Pulse Generator.

**Internal working.** The block implements a maximal-length linear feedback shift register (LFSR) of a configurable polynomial order, producing a deterministic-but-statistically-random binary sequence whose period vastly exceeds the 1024-bit sequence length configured for this run, ensuring no repetition artifacts within a single simulation execution. Because the PRBS is reseeded identically (or with a logged seed) on each sweep iteration, the sequence is reproducible for debugging and for cross-referencing against the CSV logs during verification.

**Relevant theory.** An LFSR-based PRBS approximates a Bernoulli(0.5) process over any window shorter than its period; the statistical whiteness of the sequence is what allows this single bit stream to double as the *basis-selection* random variable in the Window Selection MATLAB block without introducing a detectable bias between Z-window and X-window selection frequencies (which the SNS protocol's security proof implicitly assumes are drawn from a stable, pre-agreed probability distribution).

**Parameters and their effects.** Bit rate controls the pace at which new random bits are produced (locked to the system clock); polynomial order/seed controls sequence period and statistical properties — an insufficiently long period risks correlating the Z/X window selection pattern with the SNS Logic block's own pseudo-randomness if both ultimately derive from correlated seeds, which is a subtle implementation pitfall a rebuilder must avoid by seeding each MATLAB random-number stream independently.

**Advantages/Limitations.** PRBS generation is computationally cheap and perfectly reproducible, which aids debugging; its principal limitation is that it is *not* cryptographically secure randomness — for a simulation, this is immaterial, but a rebuilder porting this architecture toward an actual experimental control system must replace the PRBS with a certified quantum or hardware random-number generator (RNG), exactly as depicted conceptually by the "RNG" blocks in Figure 1b of Lucamarini et al.

**Signal flow / interaction with adjacent blocks.** Feeds directly into NRZ Pulse Generator, whose output electrical waveform is then read as one of the two inputs to Window Selection MATLAB.

## 4.2 User Defined Bit Sequence Generator

**Purpose.** Supplies a second, independently controllable bit sequence (labelled with an explicit bit pattern, "0101…", in the schematic) used as the payload bit source *within* Z-windows once the Window Selection stage has identified which time windows are Z-windows — i.e., this is the raw key-bit candidate stream that SNS Logic will gate according to the send/not-send decision.

**Inputs.** None (user-configured static or file-loaded bit pattern); Bit rate = Bit rate (bit/s).

**Outputs.** A single electrical bit-sequence signal to its dedicated NRZ Pulse Generator.

**Internal working.** Unlike the PRBS generator, this block plays back an explicitly specified bit pattern, which is useful during architecture verification (a rebuilder can trace a known, human-readable pattern such as alternating 0101… through the entire chain and confirm it appears correctly, and in the correct time windows, in Alice.csv) and can later be replaced with a second independent RNG stream for production-representative operation.

**Relevant theory / Parameters.** The pattern length must be commensurate with the sequence length (1024 bits) and, in production configuration, should be swapped for a properly random bit source, since in the real SNS protocol the "intent to send," not a pre-agreed bit pattern, functions as the key bit — the *value* 0 or 1 recorded in Alice's log is a direct encoding of her send/not-send decision (Step 4 of the SNS protocol, Section II of Wang et al.), not an independently chosen payload bit. In a from-scratch rebuild, this block's role should be understood as a *verification convenience*; the SNS Logic MATLAB component is what actually determines the logged Z-basis bit value, as described in Section 4.5 below.

**Advantages/Limitations.** Deterministic patterns are excellent for debugging synchronization issues but must never be mistaken for the actual security-relevant randomness of the protocol, which resides in the sending decision itself, not in this generator.

**Interaction with adjacent blocks.** Feeds its own NRZ Pulse Generator, whose output is available to the SNS Logic MATLAB block as auxiliary/verification electrical input alongside the true window-selection and sending-probability logic.

## 4.3 NRZ Pulse Generator (all instances)

**Purpose.** Converts each logical bit-sequence signal (PRBS output, user-defined sequence output, and the outputs of the Window Selection and Random Phase MATLAB blocks, each of which drives its own dedicated NRZ Pulse Generator instance in this schematic) into a rectangular non-return-to-zero electrical waveform at the system bit rate, suitable for driving optical modulators.

**Inputs.** One electrical logical bit-sequence signal.

**Outputs.** One electrical NRZ waveform, oversampled at 32 samples/bit as set globally.

**Internal working.** For each input bit, the block outputs a constant-amplitude rectangular pulse spanning the full bit period (high level for logical 1, low level for logical 0), with no intentional pulse-shaping filter applied; rise/fall time is governed by the block's internal default (typically a small fraction of the bit period, sufficient to avoid unrealistic instantaneous transitions without materially affecting the intended on/off modulation behaviour).

**Mathematical model.** The output waveform can be written as s(t) = Σ_k b_k · rect((t − kT_b)/T_b), where b_k ∈ {0, 1} is the k-th logical bit, T_b = 1/(bit rate) is the bit period, and rect(·) is the unit rectangular pulse.

**Parameters and their effects.** Amplitude levels (typically normalized 0/1 or a configured voltage swing) set the drive levels ultimately seen by the downstream modulator's electrical input port; because several of this architecture's modulators are driven in a way that must produce a *clean, high-extinction-ratio* on/off or discrete-level optical output (this is essential for accurate intensity-class discrimination in the Decoy Estimation stage), the NRZ Pulse Generator's amplitude must be calibrated against each downstream modulator's Vπ (half-wave voltage) or its equivalent normalized drive-level convention.

**Advantages/Limitations.** NRZ pulse generation is simple and low-cost to simulate but, precisely because no pulse shaping is applied, produces a spectrally broad drive waveform; this is inconsequential for the present application (the electrical waveform never itself leaves the electrical domain — it only ever drives a modulator) but would matter if this architecture were extended to also carry a classical NRZ data channel over the same optical link.

**Interaction with adjacent blocks.** Upstream: PRBS/user-defined generator, or a MATLAB component's electrical output. Downstream: the electrical drive input of an Intensity Modulator, Mach–Zehnder Modulator, or the electrical input port of a MATLAB component that needs the waveform in oversampled electrical-signal form rather than as a raw logical bit vector.

## 4.4 Window Selection MATLAB

**Purpose.** Implements Step 1 of the SNS protocol at Alice's station: for every time window i, deterministically (or pseudo-randomly, per the pre-agreed protocol probability) classify the window as a Z-window (signal window) or an X-window (decoy window), and route Alice's subsequent processing down the corresponding branch of the schematic.

**Inputs.** Input 1 (Electrical): the PRBS bit stream (one bit per time window).

**Outputs.** Output 1 (Electrical): the Z-window branch selection/gating signal. Output 2 (Electrical): the X-window branch selection/gating signal.

**Algorithm (as specified for this build).** If the PRBS bit for window i equals 1, window i is classified as a Z-window and gating output 1 is asserted for that window (permitting the Z-window branch's downstream SNS Logic and modulator chain to act on this time slot); if the PRBS bit equals 0, window i is classified as an X-window and gating output 2 is asserted instead (permitting the X-window branch's Decoy Estimation and Random Phase chain to act).

**Internal variables logged.** For every processed time window, the component appends one row to Alice's running in-memory record (later flushed to Alice_Window.csv and ultimately merged into the consolidated Alice.csv) containing: Bit Index (the window's sequential index, 0…1023), PRBS Bit (the raw 0/1 value read from the PRBS input), Selected Window (an internal numeric or enumerated code), Window Label (the human-readable string "Z" or "X"), Timestamp/Iteration (the current sweep-iteration and sample-time stamp, for cross-run bookkeeping), and Window Number (a running count of how many Z- or X-windows have been seen so far, useful for later per-basis indexing).

**Pseudocode.**
```
for each simulation sample/window i = 0 .. N-1:
    bit = read_PRBS_input(i)
    if bit == 1:
        window_label = "Z"
        assert output1 for this window; deassert output2
    else:
        window_label = "X"
        assert output2 for this window; deassert output1
    append_row_to_AliceWindowLog(
        bit_index = i,
        prbs_bit = bit,
        selected_window = window_label,
        window_number = running_count[window_label]++
    )
end for
on simulation end:
    flush_and_close(Alice_Window.csv)
```

**Flowchart (textual).**
```
[Read PRBS bit for window i]
        |
        v
   < bit == 1 ? >
     /        \
   YES          NO
    |            |
[Z-window]   [X-window]
 assert O1    assert O2
    |            |
    +----+  +----+
         |  |
   [Log row: index, bit, label, timestamp]
         |
   [Advance to window i+1]
```

**Timing.** This component must execute once per time window (i.e., once per bit period, not once per oversampled waveform sample); if the MATLAB co-simulation component is configured to be invoked once per underlying waveform sample (32 times per bit, given the global oversampling factor), the implementation must include an internal "new window" edge-detector so that exactly one CSV row is written per bit index rather than 32 duplicate or fragmentary rows — this is the single most common implementation error in this class of architecture and is called out again in Section 4.9.

**Memory usage.** The component holds, at minimum, one open file handle and a small fixed-size struct of counters; because rows are appended incrementally rather than buffered entirely in memory, peak memory usage is O(1) in sequence length, which is important for scaling this design to much larger sequence lengths than the 1024-bit value used in this particular layout.

**Interaction with OptiSystem / data synchronization.** Output 1 and Output 2 are electrical gating signals consumed by the Z-window and X-window branch's respective downstream modulator chains (Intensity Modulator on the Z-branch; Intensity Modulator plus Phase Modulator on the X-branch); the CSV log this component produces is later joined, by Bit Index, against the SNS Logic MATLAB log, the Decoy Estimation MATLAB log, and the Random Phase MATLAB log, all of which append additional columns for the *same* bit-index rows rather than creating a separate file per stage — see Section 4.8 for the consolidated Alice CSV schema.

## 4.5 SNS Logic MATLAB (Z-window branch)

**Purpose.** Implements Step 1's signal-window sub-decision and Step 4's bit-value assignment: within windows already classified as Z-windows, decide, with probability ε, whether Alice actually transmits a signal pulse ("send") or transmits nothing ("not send"), and record the corresponding Z-basis bit value.

**Inputs.** Input 1 (Electrical): the Z-window gating signal from Window Selection MATLAB (used to determine which windows this block should act upon). Input 2 (Optical): the optical carrier from the upstream CW Laser / Mach-Zehnder Modulator chain, which this block's associated Intensity Modulator will gate.

**Outputs.** Output 1 (Optical): the gated optical pulse train — a non-vacuum pulse of intensity μ′ in windows where Alice decided "send," and vacuum (extinguished output) in windows where she decided "not send."

**Algorithm.** For each Z-window, draw a Bernoulli(ε) random variable. If the draw is "send" (probability ε), assign bit value 1, command the Intensity Modulator to pass the pulse at the configured signal intensity μ′, and record Send = 1 / Not-Send = 0. If the draw is "not send" (probability 1 − ε), assign bit value 0, command the Intensity Modulator to fully extinguish the pulse (vacuum), and record Send = 0 / Not-Send = 1. This directly implements the bit-value convention specified in Step 4 of Wang et al.: "if she has decided to send out a signal pulse, she denotes a bit value 1; if she has decided not to send, she denotes a bit value 0" (Bob's convention is the mirror image, described in Chapter 5).

**Internal variables logged.** Bit Index, Window (confirmation this window is indeed "Z", cross-checked against the Window Selection log), Send (0/1), Not Send (0/1, redundant complement of Send but logged explicitly per the specification for downstream sanity-checking), Logical Decision (the underlying Bernoulli draw, useful for auditing the empirical sending probability against the configured ε), Signal State ("SIGNAL" or "VACUUM" label mirroring Send), Timestamp.

**Mathematical model.** The sending probability ε is the free parameter that, together with the signal intensity μ′, is numerically optimized (as described in Section III of Wang et al.) to maximize the final key rate at each target distance; the two-mode state produced across a Z-window, as a function of both Alice's and Bob's independent send/not-send draws, is given by Equation set following Eq. (11) of Wang et al. — vacuum ⊗ vacuum when both choose not-send, a single-photon-containing coherent state ⊗ vacuum when exactly one of them sends, and coherent state ⊗ coherent state when both send — and it is precisely this state structure that Charlie's beam-splitter measurement acts upon.

**Pseudocode.**
```
for each Z-window i (as flagged by Window Selection):
    draw = bernoulli(epsilon)
    if draw == SEND:
        bit_value = 1
        command_intensity_modulator(window=i, level=mu_prime)
    else:
        bit_value = 0
        command_intensity_modulator(window=i, level=VACUUM)
    append_row_to_AliceLog(
        bit_index = i, window = "Z",
        send = (draw==SEND), not_send = (draw!=SEND),
        logical_decision = draw, signal_state = draw,
        bit_value = bit_value
    )
end for
```

**Flowchart (textual).**
```
[Z-window i confirmed] -> [Draw Bernoulli(epsilon)]
                                  |
                        < draw == SEND ? >
                          /             \
                        YES              NO
                         |                |
              [bit=1, modulator=mu']  [bit=0, modulator=VACUUM]
                         |                |
                         +--------+-------+
                                  |
                        [Append row to Alice CSV]
```

**CSV update process.** This block appends *new columns* (Send, Not Send, Logical Decision, Signal State) onto the *same row* already created by Window Selection MATLAB for bit index i, rather than creating an independent file; achieving this in practice requires either (a) the MATLAB components sharing a common, sequentially-appended CSV via careful column-ordering discipline, with each stage writing only its own columns and leaving others blank until a final consolidation pass, or (b) each stage writing its own intermediate CSV keyed by Bit Index, with a lightweight consolidation MATLAB routine (or the offline Python stage) performing the join. This report documents both options because the schematic's naming (multiple similarly-named "MATLAB Implementation" placeholders feeding toward a single "Alice CSV" concept) is consistent with either implementation approach; Section 4.8 recommends option (b) as the more robust engineering choice.

**Parameters and their effects.** ε (sending probability): lower ε reduces the Z-basis click rate but also reduces the two-sends-collide error contribution to the Z-basis bit-error rate (since both-send and both-not-send events, which produce wrong bits, become rarer as ε shrinks — see Wang et al.'s discussion following Step 4); μ′ (signal intensity): higher μ′ increases per-pulse click probability but degrades the achievable single-photon-yield lower bound s₁ used in the decoy-state analysis; both parameters are jointly optimized in the source paper's numerical simulation, and the MATLAB component must expose both as externally configurable (ideally sweep-iteration-configurable, consistent with the "Sweep Iteration: 1/1" field visible in the layout header) rather than hard-coded.

**Advantages/Limitations.** This is the single most protocol-critical MATLAB block in the entire Alice chain: any bug that leaks the underlying Bernoulli draw's random-number-generator state, or that fails to keep Send/Not-Send decisions private until the public sifting-and-error-test stage (Steps 3–4 of the protocol) modelled downstream in Python, would constitute a genuine security flaw if this architecture were ever adapted from simulation into a real experimental controller.

**Interaction with adjacent blocks.** Consumes the Z-window gate from Window Selection MATLAB; drives the Intensity Modulator on the Z-branch; its optical output proceeds to the SMF-28 span toward Charlie (after recombination with the X-branch output, described in Section 4.7).

## 4.6 Decoy Estimation MATLAB (X-window branch)

**Purpose.** Implements the X-window portion of Step 1: for every window already classified as an X-window (decoy window), select one of several pre-agreed pulse intensity classes — full decoy intensity μ, an additional weak-decoy intensity (if a three-intensity decoy scheme is used, per the estimator in Eq. (44) of Wang et al., which explicitly requires at least two nonzero decoy intensities μ₁, μ₂ plus a vacuum reference μ₀ = 0), and vacuum — according to a configured selection probability distribution.

**Inputs.** Input 1 (Electrical): the X-window gating signal from Window Selection MATLAB. Input 2 (Optical): the optical carrier to be intensity-gated.

**Outputs.** Output 1 (Optical): the intensity-class-gated decoy pulse.

**Algorithm / Decision tree.**
```
[X-window i confirmed]
        |
[Draw intensity class from configured distribution
 {p(mu2), p(mu1), p(vacuum)}]
        |
  < which class? >
  /      |       \
mu2     mu1     vacuum
 |       |         |
[set    [set     [set
 modulator=mu2] modulator=mu1] modulator=0]
        |
[Sample mean photon number
 for logging: n_bar = selected intensity value]
        |
[Append row: bit_index, selected_intensity,
 signal(false, this is X-window), weak(bool),
 vacuum(bool), mean_photon_number, intensity_value]
```

**Internal variables logged.** Bit Index, Selected Intensity (categorical: μ₂ / μ₁ / vacuum, using the paper's μ notation for decoy intensities — note this is distinct from the Z-window's signal intensity μ′), Signal (boolean flag, false for all X-window rows, included for schema consistency with any unified downstream table), Weak (boolean, true only for the weaker of the two nonzero decoy intensities), Vacuum (boolean, true only for the μ₀ = 0 class), Mean Photon Number (the numeric μ_k value actually used for this window), Intensity Value (the corresponding physical modulator drive level / optical power setting).

**Mathematical model.** The Poisson photon-number distribution for a coherent state of mean photon number μ is p_{n|μ} = e^{−μ} μ^n / n!; the decoy-state estimator implemented downstream in Python (not in this block) uses the observed click rates S_{μ1}, S_{μ2}, and the vacuum click rate s₀, together with these p_k(μ) coefficients (explicitly, p₀(μ) = e^{−2μ}, p₁(μ) = 2μe^{−2μ}, p₂(μ) = 2μ²e^{−2μ} for the *two-mode* joint state formed when both Alice and Bob choose the same decoy intensity in a given X-window, as derived in Eqs. (6), (8), (9) of Wang et al.) to compute the lower bound on the single-photon yield s₁ (Eq. 44) and the upper bound on the single-photon phase-flip error rate ē₁^ph (Eq. 45). This MATLAB block's sole responsibility is to *generate and log* the intensity-class choice correctly and reproducibly; the estimator itself belongs entirely to the offline Python stage documented in Chapter 8.

**Pseudocode.**
```
for each X-window i (as flagged by Window Selection):
    class = sample_categorical({mu2: p2, mu1: p1, vac: p0})
    intensity_value = intensity_table[class]
    command_intensity_modulator(window=i, level=intensity_value)
    append_row_to_AliceLog(
        bit_index=i, selected_intensity=class,
        signal=false, weak=(class==mu1), vacuum=(class==vac),
        mean_photon_number=intensity_value,
        intensity_value=intensity_value
    )
end for
```

**Parameters and their effects.** Number and spacing of decoy intensities directly controls the tightness of the s₁ and ē₁^ph bounds obtained downstream; too few decoy levels (a single decoy intensity, without a vacuum reference) makes the linear system underdetermined and forces looser, more conservative bounds; too many levels dilutes the number of X-window samples available per intensity class for a fixed total sequence length, increasing statistical (finite-key) uncertainty per class even as the bound formula itself becomes tighter in the asymptotic limit. Selection probabilities {p(μ₂), p(μ₁), p(vacuum)} should sum to 1 and are, like ε and μ′, parameters that the source paper numerically optimizes per target distance.

**Advantages/Limitations.** Because this block only ever gates *decoy* windows, and because the global phase of these pulses is (unlike Z-window pulses) permitted to be publicly disclosed later, there is no security requirement to keep the intensity-class draw private in the way the Z-window Send/Not-Send draw must be kept private — this is an important asymmetry for a rebuilder to internalize, and it is the reason this MATLAB component's CSV columns can, in principle, be safely revealed to Charlie or even to Eve without compromising protocol security, exactly as Wang et al. argue at length in their Theorem (Section IV).

**Interaction with adjacent blocks.** Consumes the X-window gate; drives the Intensity Modulator on the X-branch; its output proceeds to the Phase Modulator (Section 4.7) before recombination and launch toward Charlie.

## 4.7 Random Phase MATLAB (X-window branch)

**Purpose.** Implements the global random phase δ_{Ai} required by Step 0 of the SNS protocol for every time window (in this schematic's simplified realization, applied specifically on the X-window branch where the phase value has direct, protocol-relevant consequences for the interference-based error-rate estimate; note that in the full protocol every window, Z or X, technically carries an independent random phase, but only the X-window phase is ever used in a post-selection/interference criterion).

**Inputs.** None required beyond an internal RNG seed (the block is shown with only an electrical output port in the schematic, consistent with it being a pure phase-value generator whose output drives the Phase Modulator's electrical control input).

**Outputs.** Output 1 (Electrical): a phase-control signal representing the drawn random phase value δ, scaled to the Phase Modulator's electrical-to-optical-phase transfer characteristic (typically a voltage proportional to the desired radian phase shift, calibrated against the modulator's Vπ).

**Algorithm.** For each window, draw δ uniformly from the semi-open interval [0, 2π), exactly as specified in the "Step 0" description common to both source papers ("they take random phase shifts δ_Ai, δ_Bi… taken privately by Alice and Bob"). Optionally, and consistent with the discretized phase-slice construction of Lucamarini et al. (Fig. 1c), the continuous draw may additionally be quantized into one of M phase slices Δ_k = 2πk/M for logging and later post-selection purposes; because the SNS protocol (per Wang et al.'s Note 3, Virtual Protocol 1) can operate with the phase slice made arbitrarily fine (formally λ → 0 in Eq. (1) of that paper) since no post-selection is applied in the Z-basis, this block should expose the slice count M as a configurable parameter defaulting to a fine discretization for the X-window use case, where the post-selection criterion 1 − |cos(δ_A − δ_B)| ≤ |λ| (Eq. 1 of Wang et al.) is actually applied downstream in Python.

**Internal variables logged.** Bit Index, Random Phase (the raw continuous value, radians), Phase Index (the discretized slice index k, if quantization is enabled), Phase Value (the slice-representative value 2πk/M, if applicable), Timestamp.

**Importance in SNS-TFQKD.** The random global phase is what allows a coherent-state pulse to be treated, after averaging over the phase, as a classical Poissonian mixture of photon-number states rather than as a fixed-phase coherent state with well-defined amplitude — this is the mathematical device (Eq. (5)–(9) of Wang et al.) that lets the standard decoy-state formalism, developed originally for phase-randomized BB84, apply directly to the X-basis pulses of the SNS protocol. Getting this phase draw wrong (e.g., correlating it across windows, or failing to draw it independently at Alice and at Bob) would silently invalidate the entire decoy-state estimator used downstream, without necessarily producing any visible symptom in the OptiSystem-level waveforms.

**Pseudocode.**
```
for each window i:
    delta = uniform(0, 2*pi)
    if quantization_enabled:
        k = floor(delta / (2*pi/M))
        phase_value = k * (2*pi/M)
    else:
        phase_value = delta
    command_phase_modulator(window=i, phase=phase_value)
    append_row_to_AliceLog(
        bit_index=i, random_phase=delta,
        phase_index=k, phase_value=phase_value
    )
end for
```

**Flowchart (textual).**
```
[Window i] -> [Draw delta ~ U(0, 2*pi)]
                    |
         < quantization enabled? >
             /               \
           YES                NO
            |                  |
   [k = floor(delta / slice)]  |
   [phase_value = k*slice]     [phase_value = delta]
            |                  |
            +--------+---------+
                     |
        [Drive Phase Modulator with phase_value]
                     |
        [Append row to Alice CSV]
```

**Advantages/Limitations.** Uniform phase randomization is straightforward to simulate but, in a real hardware realization, requires a fast, high-resolution phase modulator and a correspondingly fast RNG-to-drive-voltage DAC path; simulation-side, the main risk is inadvertent correlation between Alice's and Bob's independently-instantiated Random Phase MATLAB blocks if both are seeded from a shared global RNG state rather than independent seeds — a rebuilder must explicitly verify seed independence.

**Interaction with adjacent blocks.** Drives the Phase Modulator on the X-branch, positioned after the Decoy Estimation MATLAB block's Intensity Modulator in the signal-flow (intensity class is selected first, then the phase is imprinted onto whatever intensity-class pulse resulted), consistent with the schematic's left-to-right ordering: CW Laser → (Decoy Estimation) Intensity Modulator → Phase Modulator (driven by Random Phase MATLAB) → X Coupler recombination stage.

## 4.8 The Consolidated Alice CSV

**Purpose.** To provide the offline Python stage with a single, bit-index-synchronized table capturing every classical variable Alice's station produces across an entire simulation run, so that Python never needs to re-derive protocol state from raw waveforms.

**Recommended schema.**

| Column | Type | Populated by | Meaning |
|---|---|---|---|
| Iteration | int | Window Selection | Sweep-iteration index (for multi-run studies) |
| Bit Index | int | Window Selection | Global time-window index, 0…N−1 |
| PRBS Bit | 0/1 | Window Selection | Raw bit used for Z/X classification |
| Window | Z / X | Window Selection | Basis classification for this window |
| Window Label | string | Window Selection | Redundant human-readable copy of Window |
| Send | 0/1 | SNS Logic | 1 if Alice transmitted a signal pulse (Z-windows only) |
| Not Send | 0/1 | SNS Logic | Complement of Send (Z-windows only) |
| Intensity Class | μ′ / μ₂ / μ₁ / vacuum | SNS Logic / Decoy Estimation | Which intensity was actually used |
| Mean Photon Number | float | Decoy Estimation | Numeric μ value for this window (X-windows) |
| Random Phase | float (rad) | Random Phase | Raw phase draw |
| Phase Index | int | Random Phase | Quantized phase-slice index |
| Optical Power | float (dBm/W) | derived | Logged actual modulator output level, for cross-check against theoretical intensity |
| Timestamp | float | all | Simulation-time stamp of this window |

**CSV synchronization across stages.** Because Window Selection, SNS Logic, Decoy Estimation, and Random Phase are four *separate* MATLAB Component instances in the OptiSystem schematic, each executing on its own invocation schedule, the most robust implementation writes four separate intermediate files — `Alice_Window.csv`, `Alice_SNS.csv`, `Alice_Decoy.csv`, `Alice_Phase.csv` — each keyed by Bit Index, and performs a single left-join on Bit Index as the very last step of the OptiSystem run (either in a fifth, purely-bookkeeping MATLAB component that runs once at simulation end, or, more simply, in the first cell of the offline Python notebook). This avoids the substantial engineering complexity of having four independent MATLAB processes safely interleave writes to one shared file handle, at the modest cost of an extra join step that Python performs essentially instantaneously for tables of this size.

## 4.9 CW Laser, Mach–Zehnder Modulator, Intensity Modulators, Phase Modulators, Attenuators, Power Meters, Visualizers

**CW Laser (Analytical).** Purpose: provides the continuous-wave optical carrier at 193.1 THz, 0 dBm, with a configurable linewidth (relative intensity noise and phase noise parameters, left at analytical-model defaults in this build for tractability). Internal working: emits a constant-amplitude, constant-frequency complex envelope E(t) = √P₀ · e^{jφ₀}, with P₀ the configured power and φ₀ an arbitrary fixed reference phase (the *random* global phase used by the protocol is applied later, by the Phase Modulator, not by the laser itself). Interaction: feeds the Mach–Zehnder Modulator (Z-branch) or directly the Decoy Estimation Intensity Modulator (X-branch), depending on branch.

**Mach–Zehnder Modulator (Analytical).** Purpose: imprints the NRZ bit pattern from the PRBS/user-defined bit chain as an amplitude (on/off) pattern onto the optical carrier, prior to the SNS Logic stage's further intensity gating. Internal working: implements the standard dual-arm interferometric transfer function T(V) = cos²(π V / (2V_π) + φ_bias), producing high extinction ratio on/off keying when biased at quadrature and driven between 0 and V_π. Parameters: extinction ratio (must be high, >20 dB, so that "off" truly approximates vacuum rather than a dim but nonzero leakage pulse, since leakage directly inflates the effective vacuum-state click rate s₀ used in the decoy-state estimator); insertion loss (subtracted directly from the effective 0 dBm laser reference, reducing the power budget available downstream). Interaction: feeds the SNS Logic MATLAB block's associated Intensity Modulator.

**Intensity Modulator (multiple instances).** Purpose: performs the actual, MATLAB-block-commanded intensity-class gating — full extinction (vacuum) or one of the configured intensity levels (μ′ on the Z-branch; μ₂/μ₁/vacuum on the X-branch). Internal working: modelled as a linear, externally-controlled attenuator/gate whose transmittance is set, sample-by-sample, from the MATLAB Component's output rather than from a fixed OptiSystem parameter — this is the mechanism by which discrete protocol decisions (computed in MATLAB) become continuous, physically simulated optical power levels. Parameters: on/off extinction ratio (as above, critical for a clean vacuum state); response time (must be fast relative to the 100 ps bit period; an intensity modulator with response time comparable to or slower than the bit period will smear adjacent windows' intensity levels into one another, corrupting the decoy-state class boundaries). Advantages/Limitations: modelling gating as an ideal, instantaneous, MATLAB-commanded transmittance is a simplification relative to a real intensity modulator's finite bandwidth and driver-amplifier settling time; a rebuilder wanting a more realistic model should insert an explicit bandwidth-limiting filter between the MATLAB Component's electrical control output and the Intensity Modulator's control input.

**Phase Modulator.** Purpose: imprints the Random Phase MATLAB block's drawn phase value onto the X-branch optical pulse (and, in a fuller realization, would also be used at Charlie's station to apply active phase-drift compensation, as in Lucamarini et al.'s Fig. 1b, though that compensation loop is not separately modelled as a distinct block in this particular schematic beyond the phase modulators already shown at Charlie's interference stage). Internal working: applies a pure phase shift E_out(t) = E_in(t)·e^{jΔφ(t)} with Δφ(t) set by the electrical control input, ideally with negligible residual amplitude modulation (a "chirp-free" phase modulator model). Parameters: V_π (phase modulator half-wave voltage, calibrating the electrical-drive-to-phase-shift scale factor that the Random Phase MATLAB block's output must be pre-scaled against).

**Variable Optical Attenuators / Power Meters / Visualizers.** Although not each individually elaborated with its own named instance in every branch of this particular schematic (some of this functionality is folded directly into the Intensity Modulator's MATLAB-commanded transmittance, in lieu of separate discrete VOA blocks), the architectural role of a VOA, where present, is to set the *coarse* average output power of a branch (bright classical-regime reference pulses versus dim quantum-regime pulses, per Lucamarini et al.'s Fig. 1b), while Power Meters and Dual Port Oscilloscope/Binary Sequence Visualizers serve purely as engineering-verification instrumentation — they do not feed back into any protocol decision and exist solely so that a rebuilder debugging the architecture can visually confirm that the correct waveform (correct extinction ratio, correct pulse timing, correct eye pattern) appears at each stage of the chain before trusting the corresponding MATLAB CSV log.

---

# Chapter 5 — Bob Transmitter

Bob's transmitter chain is architecturally a mirror image of Alice's, instantiated as a second, independent set of blocks (Pseudo-Random Bit Sequence Generator_1, NRZ Pulse Generator_2, User Defined Bit Sequence Generator_3, NRZ Pulse Generator_7, CW Laser_1, MZ Modulator Analytical_1, Window Selection – BOB, SNS MATLAB, Intensity Modulator_2, Decoy Estimation – MATLAB [second instance], CW Laser, Phase Modulator_1, Intensity Modulator_3) feeding its own dedicated pair of SMF-28 spans toward Charlie. Every design principle, algorithm, pseudocode structure, and CSV schema documented in Chapter 4 for Alice applies identically to Bob, subject to the following protocol-specified asymmetries that a rebuilder must implement precisely rather than merely by symmetric copy-paste.

## 5.1 Bit-Value Convention Asymmetry

Per Step 4 of the SNS protocol (Wang et al.), Bob's bit-value assignment is the *logical complement* of Alice's: "if he has decided to send out a signal pulse, he denotes a bit value 0; if he has decided not to send, he denotes a bit value 1." This is not a cosmetic difference — it is what makes a "both send" or "both not-send" event in a given Z-window register as a *bit error* (Alice's bit ≠ Bob's bit is the desired outcome for a valid key bit; both agreeing is what should occur only in error) once Charlie's single-click heralding event is correlated against both users' records downstream in Python. Bob's SNS Logic MATLAB component must therefore implement:

```
if draw == SEND:
    bob_bit_value = 0     # note: opposite convention to Alice
else:
    bob_bit_value = 1
```

A rebuilder who copies Alice's SNS Logic MATLAB block verbatim for Bob's branch, without flipping this convention, will silently corrupt every downstream Z-basis error-rate calculation, because Alice's and Bob's logged bit values will systematically disagree with the sign convention the Python sifting stage expects.

## 5.2 Independent Randomness Requirement

Bob's Window Selection MATLAB, SNS Logic MATLAB, Decoy Estimation MATLAB, and Random Phase MATLAB components must each draw from random-number-generator streams that are statistically independent of Alice's corresponding streams — both because the SNS protocol's security proof explicitly treats Alice's and Bob's decisions as independently and locally generated (Step 1: "Alice (Bob) independently determines whether it is a decoy window or a signal window"), and because any simulated correlation between the two sides' random draws (e.g., from an accidentally shared global MATLAB RNG seed reused across both stations' component instances) would produce an artificially inflated or deflated Z-basis agreement rate that does not correspond to any real physical mechanism, silently distorting the simulated QBER.

## 5.3 Bob's Fibre Spans and Symmetry to Charlie

Bob's optical output, like Alice's, propagates over two independent 40 km SMF-28 spans (SMF-28 1_2 and SMF-28 1_3 in the schematic's naming) toward Charlie's station, arriving at the *same* pair of X Couplers that Alice's fields arrive at — X Coupler serving the Z-window interference path, X Coupler_1 serving the X-window interference path. Architecturally, Bob's station being placed at an equal 40 km distance from Charlie as Alice mirrors the symmetric-distance configuration used throughout both source papers' numerical simulations (equal L on each side, for a total Alice–Charlie–Bob span of 2L); a rebuilder wanting to model an asymmetric deployment (Alice closer to Charlie than Bob) must adjust the two sides' fibre-length parameters independently and should expect the achievable key rate to be governed by the *weaker* (lower-transmittance) of the two links, since Charlie's interference-based heralding event requires a successful arrival from both sides simultaneously.

## 5.4 Bob.csv

Bob's station accumulates its own consolidated CSV artifact, Bob.csv, with an identical column schema to Alice.csv (Section 4.8), populated by Bob's own instances of the Window Selection, SNS Logic, Decoy Estimation, and Random Phase MATLAB components, joined on Bit Index exactly as described for Alice. Bob.csv and Alice.csv share the same Bit Index numbering scheme (both stations process the same global sequence of time windows, 0…N−1, on the same system clock), which is what makes a straightforward Bit-Index join between the two files valid in the offline Python stage — see Chapter 7 for the full three-way (Alice/Bob/Charlie) synchronization discussion.

## 5.5 Optical and Electrical Components — Delta from Alice's Chain

All optical and electrical component types on Bob's chain (CW Laser_1, MZ Modulator Analytical_1, Intensity Modulator_2, Intensity Modulator_3, Phase Modulator_1) are configured identically to their Alice-side counterparts documented in Section 4.9 — same 193.1 THz carrier frequency, same 0 dBm reference power, same extinction-ratio and V_π considerations — with the sole intent of keeping the two transmitters symmetric so that any residual asymmetry observed in the simulated results (interference visibility, click-rate imbalance) can be attributed to the fibre-length/loss configuration or to genuine protocol statistics, rather than to an inadvertent hardware-model mismatch between the two stations. A rebuilder deliberately studying laser-frequency-offset effects (the Δν term in Eq. (4) of Lucamarini et al., governing the phase-locking requirement between Alice's and Bob's independent lasers) should detune CW Laser_1 relative to CW Laser_15 by a small, explicitly configured Δν and verify that the resulting phase-drift contribution to the X-basis error rate matches the analytical prediction of Eq. (4) before attributing any observed error to other causes.

---

# Chapter 6 — Charlie Measurement Station

## 6.1 Role of Charlie

Charlie's station is the untrusted (or semi-trusted) interference and measurement point at which Alice's and Bob's optical fields, having each propagated 40 km, are combined and jointly measured. Consistent with the measurement-device-independent character of the protocol, everything downstream of the beam splitter at Charlie's station is modelled, in this architecture, purely as a passive optical interference stage followed by click-logging detectors — no protocol decision, sifting, or key-relevant computation occurs at Charlie in this design, which is the correct architectural stance given that Charlie is never assumed trustworthy in the security analysis.

## 6.2 SMF Propagation Into Charlie's Station

The four SMF-28 spans (SMF-28 1_1, SMF-28 1, SMF-28 1_2, SMF-28 1_3 — two carrying Alice's Z- and X-branch outputs respectively, two carrying Bob's) apply standard single-mode-fibre attenuation (default OptiSystem SMF-28 attenuation coefficient, which should be explicitly verified/set to 0.2 dB/km to match the assumption used in both source papers' numerical work), chromatic dispersion, and (if enabled) polarization-mode-dispersion and nonlinear-effects models. At 40 km and the modest optical powers used in this design (well below the nonlinear threshold for a single-channel, non-WDM link), nonlinear effects are negligible and the fibre acts, to good approximation, as a linear attenuator plus a group-delay/dispersion element; the loss experienced is η_fibre = 10^{−αL/10} = 10^{−(0.2)(40)/10} ≈ 0.126, i.e., roughly 9 dB of loss per span, consistent with the channel-transmittance figure that the decoy-state and SNS key-rate formulas (Eq. 2–4 of both papers) treat as the governing loss parameter L (via η).

## 6.3 Interference at the X Coupler

**Purpose.** The X Coupler (instantiated twice — once for the Z-window path, once for the X-window path) implements the 50:50 beam-splitter interference central to both TF-QKD and SNS-TFQKD: Alice's and Bob's incoming fields are combined such that the two output ports register constructive or destructive interference depending on the *relative* phase between the two inputs, exactly as depicted in Fig. 2c of Lucamarini et al. and referenced throughout Wang et al.'s Step 2 and its associated post-selection criterion (Eq. 1).

**Internal working / mathematical model.** For two input fields with complex amplitudes E_A and E_B (each dressed with its respective random phase, encoding phase, and accumulated fibre phase noise), a standard 50:50 coupler produces output-port amplitudes E_out1 = (E_A + E_B)/√2 and E_out2 = (E_A − E_B)/√2 (up to the coupler's specific phase convention). The resulting output *intensities*, and therefore the relative click probability on each of the two downstream detectors, depend on cos(δ_A − δ_B + Δφ_encoding) — this is the physical mechanism by which "which detector clicks" carries information about whether Alice's and Bob's Z-basis send/not-send decisions were correlated in the desired sense, and, in the X-basis, about whether the phase-slice post-selection criterion (Eq. 1 of Wang et al.) is satisfied for a given pair of decoy pulses.

**Why two separate X Couplers.** The schematic's duplication of the interference stage into a Z-window instance (X Coupler, feeding SPD/SPD_1) and an X-window instance (X Coupler_1, feeding SPD_2/SPD_3) reflects the architectural decision to keep the Z-basis and X-basis optical processing paths fully separated all the way from each transmitter's Window Selection stage through to detection, rather than time-multiplexing both bases through a single shared interference stage. This duplication simplifies the MATLAB logging logic (each detector pair's click stream can be unambiguously attributed to one basis without needing an additional window-type lookup at the detection stage) at the cost of doubling the detector hardware count relative to a design that time-multiplexes a single interference stage — a legitimate engineering trade-off that a rebuilder should note when costing out a physical realization of this architecture.

## 6.4 Single-Photon Detectors (SPD, SPD_1, SPD_2, SPD_3)

**Purpose.** Convert the weak optical field at each of the two X Coupler output ports (per basis) into a discrete click/no-click electrical detection event, modelling the behaviour of a real gated avalanche photodiode or superconducting nanowire single-photon detector.

**Inputs.** One optical signal (one X-Coupler output port) per detector instance.

**Outputs.** One electrical click-event signal (and, depending on the specific OptiSystem SPD model configuration, an associated timing/jitter and dark-count-annotated output).

**Internal working.** The SPD block probabilistically registers a "click" according to a model combining (a) the instantaneous optical power at its input, converted to an expected photon count via the configured detection efficiency η_det, (b) a Poissonian dark-count process with configurable dark-count probability P_dc per gate/time-window, and (c) (if modelled) detector dead-time and afterpulsing effects that suppress or bias clicks immediately following a previous click. The click probability for a given time window, given a mean incident photon number n, is approximately P_click ≈ 1 − (1 − P_dc)·e^{−η_det·n} for a simple non-photon-number-resolving detector model, which is the standard approximation used in both source papers' analytical key-rate formulas.

**Parameters and their effects.** η_det (detection efficiency): directly scales the achievable click rate and therefore the raw key rate; both source papers use realistic values well below unity (30% in Lucamarini et al.'s "realistic" scenario, 80% in Wang et al.'s SNS numerical simulation, reflecting different detector technology assumptions — a rebuilder should choose a value consistent with the detector technology (InGaAs APD versus superconducting nanowire) being modelled). P_dc (dark-count probability): directly sets a floor under which real signal cannot be distinguished from detector noise, becoming the dominant error contribution at long distances where the true signal click rate falls toward or below the dark-count rate — this is, physically, what ultimately limits the maximum secure distance in both papers' numerical results. Dead-time: limits the maximum sustainable click rate and, at very high bit rates, can introduce click-probability correlations between adjacent time windows that a rigorous finite-key security analysis must account for (though, consistent with the asymptotic scope of both source papers' headline results, this architecture's default configuration treats each window's detection event as statistically independent).

**Advantages/Limitations.** Modelling detection purely as a probabilistic click/no-click event, without any photon-number-resolving capability, is standard and matches both papers' theoretical treatment (which explicitly assumes non-photon-number-resolving detectors and instead relies on the decoy-state method to statistically infer single-photon contributions); a rebuilder wanting to model a photon-number-resolving detector technology would need to extend the SPD block's output to a small integer photon-count value rather than a binary click flag, and would need to correspondingly extend every downstream MATLAB logging component and the offline Python analysis to carry that additional information.

## 6.5 SPD Detection Logging

**Purpose.** Capture, for every time window and every one of the four detectors, the full detection-event record needed for downstream sifting, effective-event determination, and error-rate estimation.

**Required log fields per detector, per window.** Bit Index, Detector Name (D0/D1 for the Z-window pair, D2/D3 for the X-window pair, or an equivalent labelling scheme distinguishing all four SPD instances), Click (boolean), No Click (boolean, complement of Click), Arrival Time (the detection timestamp within the gate/window, useful for jitter analysis and for distinguishing genuine signal-correlated clicks from dark counts in a more detailed timing-resolved analysis), Detection Probability (the instantaneous modelled click probability for this window, useful for verification against the analytical model above), Dark Count (boolean flag indicating whether the click, if any, is attributable — in the simulation's internal random-process bookkeeping — to the dark-count process rather than to genuine incident light; note that in a real experimental system this flag would not be knowable and is included here purely as a simulation-debugging aid, not as something the offline Python analysis is permitted to use, since a real analysis pipeline must estimate dark-count contamination statistically rather than by direct labelling), Detection Event (a consolidated categorical summary — "no event," "single click D_x," "double click" (i.e., both detectors of a pair clicking in the same window, which per the SNS protocol's effective-event definition invalidates that window rather than heralding a valid event)).

**Pseudocode for the logging component.**
```
for each time window i:
    for each detector d in [D0, D1, D2, D3]:
        click, is_dark_count, arrival_time = query_SPD(d, i)
        append_row_to_DetectorLog(
            bit_index=i, detector_name=d,
            click=click, no_click=!click,
            arrival_time=arrival_time,
            dark_count=is_dark_count
        )
    # after collecting both detectors of a pair, classify the window:
    pair = (basis==Z) ? [D0,D1] : [D2,D3]
    n_clicks = count_clicks(pair, window=i)
    if n_clicks == 1:
        event = "single-click:" + clicking_detector_name
    elif n_clicks == 0:
        event = "no-click"
    else:
        event = "double-click (invalid)"
    append_event_summary(bit_index=i, basis=basis, event=event)
end for
```

**Detector timing.** Because this architecture operates its single-photon detectors synchronously, gated once per time window (100 ps period at the configured 10 Gbit/s bit rate), the SPD logging component must, exactly like the transmitter-side MATLAB components, guard against writing multiple log rows per underlying oversampled waveform sample; the same "new window" edge-detection discipline described in Section 4.4 applies here without modification.

**Flowchart (textual, per detector pair).**
```
[Window i: query D_a, D_b]
        |
[Count clicks n = clicks(D_a) + clicks(D_b)]
        |
   < n == 1 ? >
    /       \
  YES         NO
   |           |
[valid       < n == 0 ? >
 effective     /     \
 event]      YES      NO
   |          |        |
   |     [no event] [double-click,
   |                  discard window]
   +----------+---------+
              |
     [Append event classification to Detector CSV]
```

## 6.6 Data Recovery Blocks

**Purpose.** Reconstruct a clean, synchronized logical bit/click sequence from each detector's raw electrical output, correcting for the detector model's internal timing offsets and providing a stable digital signal for the downstream Dual Port Binary Sequence Visualizer and for the CSV logging MATLAB component's window-aligned sampling. Two Data Recovery instances per basis (Data Recovery / Data Recovery_1 for the Z-window pair; Data Recovery_2 / Data Recovery_3 for the X-window pair) mirror the two detectors per basis.

**Internal working.** Implements clock-and-data recovery against the known system bit-clock (rather than an independently recovered clock, since this simulation, unlike a real experimental system, has direct access to the transmitter-side clock reference), aligning each detector's asynchronous click events onto the discrete window-index grid used throughout the rest of the architecture.

**Interaction with adjacent blocks.** Feeds the Dual Port Binary Sequence Visualizer (engineering verification only) and the CSV FINAL [Z/X]-WINDOW MATLAB logging component, which is the terminal stage responsible for producing the consolidated, synchronized detector-side CSV described in Section 6.5 and Chapter 7.

## 6.7 Visualizers

**Dual Port Oscilloscope Visualizer (per detector pair).** Displays the raw electrical waveform from each detector pair for visual inspection of pulse shape, timing, and gross signal-to-noise behaviour; purely diagnostic, contributing nothing to the CSV logs or downstream analysis.

**Dual Port Binary Sequence Visualizer (per basis).** Displays the recovered digital click/no-click sequence from each Data Recovery pair side-by-side, allowing a rebuilder to visually confirm, during architecture verification, that click events are landing on the expected window boundaries and that the double-click / single-click / no-click classification implemented in the CSV logging MATLAB component (Section 6.5) matches what is visually apparent in the recovered bit pattern.

## 6.8 CSV Final Z-Window / X-Window MATLAB (Terminal Logging Stage)

**Purpose.** The terminal MATLAB Component at Charlie's station (CSV FINAL Z-WINDOW – MATLAB, and its analogous X-window instance) consolidates the Data-Recovery-aligned click information from both detectors of a basis pair into the final, per-basis detector CSV, ready for the offline Python join against Alice.csv and Bob.csv.

**Inputs.** The recovered digital sequences from both Data Recovery blocks of the relevant basis pair.

**Outputs.** A written CSV file (conceptually Detector_Z.csv and Detector_X.csv, or a single unified Detector.csv distinguished by a Basis column) containing, per window, the full schema described in Section 6.5.

**Algorithm.** Directly implements the per-window event-classification pseudocode of Section 6.5, plus final-flush/file-close logic triggered on simulation-end (mirroring the transmitter-side MATLAB components' end-of-run flush behaviour described in Section 4.4).

---

# Chapter 7 — CSV Synchronization

## 7.1 The Three (or More) Files

By the end of one OptiSystem simulation run, the architecture has produced, at minimum: **Alice.csv** (or its four constituent intermediate files, per the recommendation in Section 4.8), **Bob.csv** (likewise), and **Detector.csv** (or its Z/X-split constituent files, per Section 6.8). Each file is indexed by the same global Bit Index / time-window numbering scheme, 0 through N−1 (N = 1024 in this layout's default configuration), because every MATLAB component in the entire architecture, at both transmitters and at Charlie, is driven off the same underlying OptiSystem sample clock (sample rate 3.2 × 10¹¹ Hz, 32 samples per bit) and therefore shares a common notion of "window i."

## 7.2 Timestamp Synchronization

Although Bit Index is the primary join key, every log row also carries a Timestamp field (the simulation-time value, in seconds, at which the corresponding window was processed). Because the physical fibre spans introduce a real propagation delay (≈196 μs one-way at 40 km, using a group velocity of ≈2.04 × 10⁸ m/s typical of SMF-28), a signal that Alice's SNS Logic MATLAB block logs at bit index i and simulation-time t does not physically arrive at Charlie's X Coupler until simulation-time t + τ_fibre; if a rebuilder's offline analysis performs any timestamp-based correlation (rather than the recommended, and much simpler, Bit-Index-based correlation), the fibre propagation delay must be explicitly added when comparing Alice's/Bob's transmit timestamps against Charlie's detection timestamps. The Bit-Index join sidesteps this issue entirely, because OptiSystem's fibre block itself correctly delays the *waveform*, meaning that by the time Charlie's MATLAB components observe "window i," they are observing the physically correct, delay-compensated instance of that window — the Bit Index label travels with the logical window, not with wall-clock simulation time, provided every MATLAB component in the chain is written to key its logs off a monotonically-incrementing logical window counter rather than off raw simulation time.

## 7.3 Index Synchronization and Packet Synchronization

A critical implementation discipline, referenced repeatedly in Chapters 4 through 6, is that every MATLAB Component in this architecture must write **exactly one CSV row per logical time window**, never one row per underlying oversampled waveform sample (32 samples per bit in this configuration) and never more than one row per window even when a component happens to be invoked multiple times within a single window's 32-sample span. The recommended implementation pattern is a persistent internal counter, incremented only when a "new window" edge is detected (e.g., by comparing the current sample's window-aligned time-stamp against the last-processed window index, and only appending a row and advancing internal protocol-decision state when the two differ), with all 32 within-window sample invocations after the first being no-ops from the standpoint of CSV writing (though they may still be necessary from the standpoint of holding a modulator's control signal constant across the whole window, which is a separate, waveform-level concern distinct from the once-per-window logging concern).

## 7.4 Data Integrity and Missing Samples

Because Charlie's detection process is inherently probabilistic (a given window may see zero clicks at all — the overwhelming majority case at any non-trivial channel loss, since η_fibre ≈ 0.126 per 40 km span combined with η_det well below unity yields a joint single-photon arrival-and-detection probability far below 1 per window), the Detector CSV will necessarily contain many more "no-click" rows than "click" rows, and it is essential that **every** window still receives a logged row (with Click = 0 for both detectors) rather than only logging rows for windows in which a click actually occurred. Omitting no-click rows would silently corrupt the offline Python stage's ability to compute correct gain values (Q = clicks / total windows sent, per basis and intensity class), since the denominator of that ratio depends on knowing the true total window count, not just the count of interesting (clicking) windows.

## 7.5 Alignment Verification

A recommended verification step, to be performed immediately after each OptiSystem run and before handing the CSV files to the main offline Python analysis, is a lightweight alignment check: confirm that Alice.csv, Bob.csv, and Detector.csv each contain exactly N = 1024 rows (or N × number-of-detectors rows for the long-format Detector CSV), confirm that the Bit Index column in each file is a contiguous, gap-free sequence 0…N−1 (any gap indicates a dropped window, most commonly caused by an off-by-one error in the "new window" edge-detection logic described in Section 7.3, or by a MATLAB component crashing partway through the run without properly flushing its file handle), and confirm that the Window / Window Label columns in Alice.csv and Bob.csv agree window-by-window (since Window Selection is driven by each station's own independent PRBS, Alice's and Bob's Z/X classification for a given window will *not*, in general, match — both source papers explicitly account for "mismatching windows" in which one party commits to a signal window while the other commits to a decoy window, per Wang et al.'s discussion of the complete SNS protocol in Section V-E — but the alignment check should confirm that this mismatch rate is statistically consistent with the configured Z/X selection probabilities, rather than symptomatic of a synchronization bug).

---

# Chapter 8 — Final Processing (External Python Environment)

## 8.1 Statement of Scope

No final key generation, no QBER computation, no secret-key-rate computation, no parameter estimation, and no finite-key analysis is performed inside OptiSystem. This is a firm architectural boundary, not a simplification of convenience: OptiSystem's role in this design ends the moment the three (or more) synchronized CSV files described in Chapter 7 have been written to disk. Everything from that point forward — sifting, matching, detector correlation, gain calculation, QBER calculation, error-correction-cost accounting, decoy-state parameter estimation, finite-key statistical bounding, secure-key-rate computation, and the production of publication-style graphs — is the responsibility of an external Python environment operating entirely offline from the physical-layer simulation.

## 8.2 CSV Loading and Synchronization

The Python pipeline begins by loading Alice.csv, Bob.csv, and Detector.csv (or their constituent intermediate files, performing the Bit-Index joins recommended in Sections 4.8 and 6.8 if the OptiSystem-side consolidation was deferred to this stage), verifying row counts and index contiguity per the checks described in Section 7.5, and producing a single master DataFrame indexed by Bit Index with columns drawn from all three sources: Alice's Window/Send/Intensity/Phase columns, Bob's equivalent columns, and Charlie's per-detector Click/Arrival-Time columns for the relevant basis.

## 8.3 Sifting and Matching

Sifting discards every window in which Alice's and Bob's window-type commitments did not match (Alice chose Z while Bob chose X, or vice versa) — these are the "mismatching windows" of Wang et al.'s Virtual Protocol 5 — retaining only genuine Z-windows (both committed to signal) and genuine X-windows (both committed to decoy, using the same intensity class, per the X-window definition in Section II of the SNS paper). Within the retained X-windows, a further post-selection is applied per Eq. (1) of Wang et al.: only those windows where the logged random-phase values satisfy 1 − |cos(δ_A − δ_B)| ≤ |λ| are kept as valid X-basis samples for interference-based error-rate estimation; all others are discarded.

## 8.4 Detector Correlation and Effective-Event Determination

For each retained window (Z or X, post-sifting), the pipeline consults Detector.csv for the corresponding basis pair and classifies the window as an *effective event* if and only if exactly one of the two detectors clicked (Section II, Step 2 of Wang et al.); windows with zero clicks are discarded as non-events, and windows with both detectors clicking are discarded as invalid (double-click) events, consistent with the effective-event definition used throughout the security proof (Virtual Protocol 1's Definition of effective event, Wang et al.).

## 8.5 Gain Calculation

For each intensity class k (including the Z-basis "signal" class and the X-basis decoy classes μ₁, μ₂, and the vacuum reference μ₀ = 0), the pipeline computes the observed gain (yield) as Q_k = (number of effective events at intensity class k) / (number of windows sent at intensity class k), directly analogous to the S_μ notation used throughout Wang et al. (e.g., S_{μ1}, S_{μ2}, s₀ in Eq. 44).

## 8.6 QBER Calculation

**Z-basis QBER.** For every effective Z-window event, the pipeline compares Alice's logged bit value against Bob's logged bit value (recalling the complementary bit-value convention of Section 5.1: Alice's "send" → 1, Bob's "send" → 0), counts a bit error whenever the two values are *equal* rather than complementary (i.e., whenever both sent or both did not send, per the discussion following Step 4 of the SNS protocol), and reports E_Z as the ratio of such errors to the total number of effective Z-window events (after removing a randomly selected error-test subsample per Step 4, whose bit values are consumed purely for the QBER estimate and discarded from the key-generating pool).

**X-basis error rate (misalignment / phase-flip proxy).** For every effective X-window event (post phase-slice selection, Section 8.3), the pipeline classifies a right versus a wrong X-bit according to the sign of cos(δ_A − δ_B) versus which detector clicked (Section II, "Note" following Step 5, Wang et al.: a right X-bit is left-detector-clicking for positive cos(δ_A − δ_B), or right-detector-clicking for negative cos(δ_A − δ_B); a wrong X-bit is the opposite pairing), and reports E^X_μ per intensity class.

## 8.7 Error Correction Preparation

The pipeline computes the classical information-reconciliation cost using the binary entropy function h(x) = −x log₂ x − (1−x) log₂(1−x) scaled by the configured error-correction efficiency factor f (f = 1.1 to 1.15 in the values used by both source papers), producing the f·H(E_Z) term that appears directly in the final key-length formula (Eq. 3, Wang et al.).

## 8.8 Parameter Estimation (Decoy-State Analysis)

Using the per-intensity-class gains and X-basis error rates from Sections 8.5–8.6, the pipeline evaluates the decoy-state single-photon yield lower bound s₁ (Eq. 44 of Wang et al.) and the single-photon phase-flip-error-rate upper bound ē₁^ph (Eq. 45), from which the number of untagged (single-photon-caused) Z-basis bits n₁ is derived via n₁ = 2ε(1−ε)μ′e^{−μ′}s₁ × (total Z-windows), consistent with the per-window key-rate formula (Eq. 4, Wang et al.).

## 8.9 Finite-Key Statistics

For a from-scratch rebuild intended to produce results directly comparable to the asymptotic curves in Figures 1–2 of Wang et al., the pipeline may initially implement only the asymptotic formulas reproduced above (as both source papers' headline numerical results are explicitly asymptotic — "we have only considered the asymptotic result," Section III of Wang et al.); a rebuilder wishing to extend the analysis to finite-key security bounds (accounting for the statistical uncertainty inherent in estimating s₁, ē₁^ph, and E_Z from a finite sample of 1024 (or more, per Section 3.1.2) simulated windows) should apply standard Chernoff-bound or Hoeffding-inequality-based confidence-interval widening to each estimated quantity before substituting into the key-length formula, at the cost of a reduced but rigorously secure key-length figure relative to the asymptotic estimate.

## 8.10 Secret Key Rate Computation

The final number of secure key bits is computed via the SNS protocol's key-length formula:

N_f = n₁ − n₁·H(ē₁^ph) − n_t·f·H(E_Z)

where n_t is the total number of retained (post-sifting, post-error-test-removal) effective Z-window events, matching Eq. (3) of Wang et al. Dividing N_f by the total simulated time (or total number of pulses sent) yields the secret key rate in bits per second (or bits per pulse), directly comparable to the log-scale key-rate-versus-distance curves of Fig. 1 in Wang et al. and Fig. 1a of Lucamarini et al.

## 8.11 Visualization and Final Graphs

The Python pipeline's final stage reproduces, for the rebuilder's own simulated parameter set, the standard SNS-TFQKD result plots: secret key rate (log scale) versus Alice–Bob distance, for a family of misalignment-error or detector-parameter sweeps (mirroring Fig. 1 of Wang et al.), and secret key rate versus a single swept parameter (e.g., misalignment error at fixed distance, mirroring Fig. 2 of Wang et al.), using the "Sweep Iteration" mechanism visible in the OptiSystem layout header to drive a batch of simulation runs across the desired parameter grid, each producing its own CSV triple to be ingested and reduced to one key-rate data point by this offline pipeline.

## 8.12 Why This Architecture (Recap)

OptiSystem, in this design, acts strictly as the physical-layer simulator: it is responsible for faithfully reproducing the classical and quantum-optical physics of laser emission, modulation, fibre propagation, interference, and single-photon detection, parameterized to match realistic (or intentionally idealized, for comparison against theoretical bounds) hardware. Python, receiving only the synchronized CSV artifacts this simulation produces, performs the entirety of the protocol-layer cryptographic analysis. This separation is what allows the same physical-layer simulation run to be re-analyzed under multiple different protocol-layer assumptions (asymptotic versus finite-key, different error-correction efficiency factors, different post-selection phase-slice widths) without ever needing to re-run the computationally expensive OptiSystem simulation — and, just as importantly, it is what keeps the physical-layer model honest, in the sense that no block inside the OptiSystem schematic is ever in a position to "cheat" by using information (such as Alice's private Send/Not-Send decision) that a real Charlie, or a real Eve, would not actually have access to during the quantum transmission phase of the protocol.

---

# Appendix A — Full Component Table

| Component | Purpose | Inputs | Outputs | Signal Type | Key Parameters | Engineering Notes |
|---|---|---|---|---|---|---|
| Pseudo-Random Bit Sequence Generator | Master Z/X selection randomness | — | Electrical bit stream | Electrical | Bit rate | Reseed independently per station |
| User Defined Bit Sequence Generator | Verification / payload bit source | — | Electrical bit stream | Electrical | Bit rate, pattern | Debug aid; real key bit comes from SNS Logic |
| NRZ Pulse Generator | Bit-to-waveform conversion | Electrical bits | Electrical NRZ waveform | Electrical | Amplitude, samples/bit | No pulse shaping filter applied |
| CW Laser | Optical carrier | — | Optical CW field | Optical | Frequency 193.1 THz, Power 0 dBm | Pre-attenuation reference power |
| MZ Modulator Analytical | NRZ amplitude modulation | Optical + electrical | Optical NRZ | Optical | Extinction ratio, V_π | High ER required for clean vacuum |
| Window Selection MATLAB | Z/X classification | Electrical PRBS | 2× electrical gate | Electrical | — | Must log once per window, not per sample |
| SNS Logic MATLAB | Send/not-send decision | Electrical gate, optical | Optical gated pulse | Mixed | ε, μ′ | Alice bit=1 on send; Bob bit=0 on send |
| Decoy Estimation MATLAB | Intensity class selection | Electrical gate, optical | Optical gated pulse | Mixed | μ₁, μ₂, probabilities | Never needs to stay private |
| Random Phase MATLAB | Global phase draw | — | Electrical phase control | Electrical | M (phase slices) | Independent seed per station |
| Intensity Modulator | Physical intensity gating | Optical, electrical control | Optical | Mixed | Extinction ratio, response time | Commanded sample-by-sample by MATLAB |
| Phase Modulator | Physical phase imprinting | Optical, electrical control | Optical | Mixed | V_π | Chirp-free assumed |
| SMF-28 Fibre | Channel propagation | Optical | Optical | Optical | Length 40 km, α = 0.2 dB/km | Adjust for long-haul studies |
| X Coupler | Beam-splitter interference | 2× optical | 2× optical | Optical | Split ratio 50:50 | Duplicated per basis |
| SPD | Single-photon detection | Optical | Electrical click | Mixed | η_det, P_dc | Non-photon-number-resolving |
| Data Recovery | Click stream alignment | Electrical | Electrical | Electrical | — | Aligns to system clock, not recovered clock |
| CSV FINAL [Z/X]-WINDOW MATLAB | Detector-side logging | Electrical (both detectors) | CSV file | — | — | Must log 0-click windows too |
| Dual Port Oscilloscope Visualizer | Diagnostic display | Electrical/optical | Display | — | — | No downstream effect |
| Dual Port Binary Sequence Visualizer | Diagnostic display | Electrical | Display | — | — | No downstream effect |

# Appendix B — Consolidated Mathematical Model

**Laser output (ideal CW):** E(t) = √P₀ · e^{jφ₀}

**Coherent state / Poisson photon statistics:** p_{n|μ} = e^{−μ}μⁿ/n!

**Fibre attenuation:** η_fibre(L) = 10^{−αL/10}, α ≈ 0.2 dB/km

**MZM transfer function:** T(V) = cos²(πV/(2V_π) + φ_bias)

**Beam-splitter outputs:** E_out1 = (E_A + E_B)/√2 , E_out2 = (E_A − E_B)/√2

**Detector click probability:** P_click ≈ 1 − (1 − P_dc)·e^{−η_det·n}

**Two-mode Z-window photon-number decomposition (Wang et al., Eq. 5–9):**
p₀(μ) = e^{−2μ}, p₁(μ) = 2μe^{−2μ}, p₂(μ) = 2μ²e^{−2μ}

**Decoy-state single-photon yield lower bound (Wang et al., Eq. 44):**
s₁ ≥ [p₂(μ₂)(S_{μ1} − p₀(μ1)s₀) − p₂(μ1)(S_{μ2} − p₀(μ2)s₀)] / [p₂(μ2)p₁(μ1) − p₂(μ1)p₁(μ2)]

**Phase-flip error-rate upper bound (Wang et al., Eq. 45):**
ē₁^ph = [S_{μ1}E^X_{μ1} − e^{−2μ1}s₀/2] / [2μ1 e^{−2μ1} s₁]

**Final key length (Wang et al., Eq. 3):**
N_f = n₁ − n₁H(ē₁^ph) − n_tf·H(E_Z)

**Per-window key rate (Wang et al., Eq. 4):**
R = 2ε(1−ε)μ′e^{−μ′}s₁[1 − H(ē₁^ph)] − S_Zf·H(E_Z)

**TF-QKD asymptotic key rate (Lucamarini et al., Eq. 3):**
R_TF−QKD^(¬ρ)(μ,L) = (d/M)[R_QKD(μ, L/2)]_⊕E_M

**Intrinsic QBER from phase-slice sifting (Lucamarini et al., Eq. 1):**
E_M = (1/2) − sin(2π/M)/(4π/M)

**Differential phase-fluctuation model (Lucamarini et al., Eq. 4):**
δ_ba = (2π/s)(ΔνL + νΔL)

**Binary entropy function:** H(x) = −x·log₂x − (1−x)·log₂(1−x)

# Appendix C — Referenced Source Papers

1. M. Lucamarini, Z. L. Yuan, J. F. Dynes, A. J. Shields, "Overcoming the rate–distance limit of quantum key distribution without quantum repeaters," *Nature* 557, 400–403 (2018).
2. X.-B. Wang, Z.-W. Yu, X.-L. Hu, "Sending or not sending: Twin-field quantum key distribution with large misalignment error," arXiv:1805.09222v9 [quant-ph] (2018).
3. SNS-TFQKD OptiSystem architecture schematic, "Layout 1," Sweep Iteration 1/1, dated 17 July 2026 (source layout documented in this report).

---

*End of report.*