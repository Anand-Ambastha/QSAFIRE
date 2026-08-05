# SNS-TF-QKD - Summer Internship, Scientific Analysis Group, DRDO

**Anand Kumar** · B.Tech, Electronics and Communication Engineering, Bharati
Vidyapeeth's College of Engineering · Summer Internship, June–August 2026,
Scientific Analysis Group (SAG), DRDO, Metcalfe House, Delhi · Supervisor:
Rajesh Kumar, Scientist 'E', SAG.

Full report: [`Summer_Internship_Report_DRDO.pdf`](./Summer_Internship_Report_DRDO.pdf)
("*QSAFire: Design and Performance Evaluation of Satellite-Assisted
SNS-TF-QKD over the Delhi-Mumbai Quantum Communication Link*").

## Three independent tasks

The internship comprised three technical tasks, all concerned with
Sending-or-Not-Sending Twin-Field Quantum Key Distribution (SNS-TF-QKD), but
carried out as separate pieces of work using separate tools — not stages of
one pipeline. No output of Task 1 or Task 2 feeds into Task 3, and this
repository keeps that boundary visible rather than merging them:

```
QSAFire-Impact-DRDO/
├── osd/                        # Task 1 + Task 2 (one path: OptiSystem -> Python pipeline)
│   ├── optisystem_prototype/     #   Task 1: OptiSystem SNS-TF-QKD prototype + MATLAB co-sim
│   └── preprocessing_pipeline/   #   Task 2: Python CSV post-processing pipeline
├── qsafire/                    # Task 3: independent satellite-link simulation framework
└── Summer_Internship_Report_DRDO.pdf
```

| | Task 1 | Task 2 | Task 3 |
|---|---|---|---|
| Folder | `osd/optisystem_prototype/` | `osd/preprocessing_pipeline/` | `qsafire/` |
| What it is | Physical-layer SNS-TF-QKD link in OptiSystem, protocol logic via embedded MATLAB co-simulation blocks | Python pipeline that ingests Task 1's CSV exports and runs the Wang/Yu/Hu decoy-state + key-rate analysis | Independent Python framework evaluating SNS-TF-QKD over a satellite-to-ground downlink |
| Channel | 40 km SMF-28 fibre spans, laboratory scale | (consumes Task 1's output only) | Satellite-to-ground free-space, 234–1,150 km |
| Report chapter | Ch. 2 | Ch. 3 | Ch. 4–5 |

**Task 1 and Task 2** form one path — OptiSystem output consumed by the
Python pipeline — and live together under `osd/` per the report's own
naming. **Task 3, QSAFire,** shares no code, data, or CSV schema with either
and is fully self-contained in `qsafire/`.

### Task 1 — OptiSystem SNS-TF-QKD Prototype (`osd/optisystem_prototype/`)

A three-station (Alice / Bob / Charlie) proof-of-concept link built in
OptiSystem: fibre spans, CW lasers, Mach-Zehnder and intensity modulators,
and single-photon detectors at the physical layer, with the SNS protocol's
decision logic (window selection, the send/not-send decision, decoy-intensity
selection, phase randomization) implemented through embedded MATLAB
co-simulation blocks. Laboratory scale — 40 km spans, 10 Gbit/s, a 1024-bit
sequence per run — with no satellite or orbital modeling of any kind. See
report Chapter 2.

### Task 2 — Python Preprocessing Pipeline (`osd/preprocessing_pipeline/`)

Reads the CSV logs Task 1 exports, synchronizes them by window index, and
runs the classical post-processing the SNS security analysis requires:
basis separation, X-basis phase post-selection, decoy-state yield/phase-error
estimation (Eqs. 44–45), Z-basis gain/QBER, and the secret-key-rate formula
(Eq. 4) — with a Wilson confidence interval on every finite-count estimate.
Operates only on OptiSystem-generated data. See report Chapter 3.

### Task 3 — QSAFire (`qsafire/`)

*Quantum Satellite Analysis Framework for Integrated Research and
Evaluation.* An independent, from-first-principles simulation of an
SNS-TF-QKD uplink to a sun-synchronous satellite: orbital geometry (inclination
solved from J2 perturbation theory, validated against Landsat-8/9 and
Sentinel-2), atmospheric transmittance, slant-path Rytov-variance turbulence
(including a corrected spherical-wave/uplink formula), pointing loss, real
hardware/background-light losses, weather-gated annual availability, joint
signal/decoy optimization, and a fully asymmetric per-arm secure-key-rate
calculation across six Indian ground sites and four candidate links. The
primary Delhi–Mumbai link returns a simulated secure key rate of
0.0947 kbit/s over 1,149.4 km — roughly two orders of magnitude below
Micius's flagship figures but consistent with real, non-ideal portable
ground-station demonstrations in the literature — and identifies night-sky
background light, not turbulence, as the binding constraint. See report
Chapter 4–5, and `qsafire/README.md`, `qsafire/docs/api.md`, and
`qsafire/docs/theory.md` for full documentation of this component.

## Adding the Task 1 / Task 2 files

This repository was assembled with `qsafire/` as a complete, executable
package (see its own README for setup and usage) and the internship report
as the authoritative record of Tasks 1 and 2. The Task 1 OptiSystem
project file, its MATLAB co-simulation sources, and the Task 2 pipeline
notebook/script and CSV inputs were not available when this structure was
put together — `osd/README.md` lists exactly where each belongs and what
each is expected to contain, drawn directly from the report.

## Where to start

- **Reading the work end to end:** `Summer_Internship_Report_DRDO.pdf`.
- **Running the satellite-link simulation:** `qsafire/README.md`.
- **Understanding QSAFire's physics:** `qsafire/docs/theory.md`.
- **QSAFire's function-level API:** `qsafire/docs/api.md`.
- **Task 1 / Task 2 file layout and what belongs where:** `osd/README.md`.
