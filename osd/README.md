# OSD — OptiSystem SNS-TF-QKD Prototype & Preprocessing Pipeline

This folder holds **Task 1** and **Task 2** of the internship: a physical-layer
SNS-TF-QKD prototype built in OptiSystem (with protocol logic implemented as
MATLAB co-simulation blocks), and the Python pipeline that post-processes the
CSV logs that prototype exports. Together they form one path — OptiSystem
output consumed by the Python pipeline — and are entirely separate from
`../qsafire`, which shares no code, data, or CSV format with either of them.

Full documentation for both tasks is in
`../Summer_Internship_Report_DRDO.pdf`, Chapters 2 and 3. This README is a
navigation aid, not a replacement for that report.

## `optisystem_prototype/` — Task 1

A three-station (Alice / Bob / Charlie) SNS-TF-QKD link built in OptiSystem,
using single-mode fibre (SMF-28) spans between stations and a 50:50 coupler
+ single-photon-detector pair at Charlie for each basis. The protocol-level
decision logic — per-window Z/X basis classification, the Z-window
send/not-send decision, X-window decoy-intensity selection, and X-window
phase randomization — has no native OptiSystem block, so it is implemented
as MATLAB co-simulation components acting as programmable taps in the
signal-flow graph (see report Table 2.1 for the role of each block).

**Expected contents** (add the actual project files here):

```
optisystem_prototype/
├── <project_name>.osd          # the OptiSystem schematic/project file
├── matlab_cosim/                # MATLAB co-simulation block source
│   ├── window_selection.m
│   ├── sns_logic_alice.m
│   ├── sns_logic_bob.m           # NOTE: bit-value convention is the
│   │                              #   logical complement of Alice's —
│   │                              #   see report Section 2.2
│   ├── decoy_estimation.m
│   ├── random_phase.m
│   └── csv_logging.m
└── outputs/                     # exported CSVs from a simulation run
    ├── alice_signal_record.csv
    ├── bob_signal_record.csv
    ├── charlie_record_x.csv
    ├── charlie_record_z.csv
    └── charlie_xbasis_log.csv
```

Key simulation parameters (report Table 2.2): 10 Gbit/s bit rate, 1024-bit
sequence per run, 32 samples/bit, 193.1 THz (≈1552.52 nm) CW laser, 0 dBm
launch power, 4× 40 km SMF-28 spans (≈9 dB loss/span at 0.2 dB/km).

**The `.osd` project file, the MATLAB `.m` sources, and the raw CSV outputs
were not included when this folder structure was set up — add them here
directly from the OptiSystem workstation.**

## `preprocessing_pipeline/` — Task 2

A Python pipeline that reads the five CSV files exported by the Task 1
prototype (`alice_signal_record.csv`, `bob_signal_record.csv`,
`charlie_record_x.csv`, `charlie_record_z.csv`, `charlie_xbasis_log.csv`),
synchronizes them by `Window_Index`, and carries out the classical
post-processing Wang, Yu & Hu's SNS security analysis requires:

1. **Synchronization** — row-count / duplicate-index / missing-column checks
   across all five files.
2. **Basis separation** — discard windows where Alice's and Bob's
   independently-chosen Z/X basis labels don't match.
3. **X-basis phase post-selection** — decoy-intensity matching, then the
   phase-acceptance rule `1 - |cos(δ_A - δ_B)| ≤ |λ|` (Eq. (1) of Wang et
   al.), applied via the discretized phase-slice convention the OptiSystem
   model uses.
4. **Decoy-state estimation** — the two-decoy-intensity yield lower bound
   (Eq. 44) and phase-error upper bound (Eq. 45).
5. **Z-basis gain and QBER** — split into the four SNS categories (Alice
   sends / Bob sends / both send / neither sends); only the first two carry
   genuine key bits.
6. **Secret key rate** — Wang et al.'s Eq. (4), with every constituent term
   reported and a 95% Wilson score confidence interval on every rate
   estimated from a finite count.

See report Section 3.9 for an honest read of the one CSV set this pipeline
was run against: it executed end to end without error and returned a
negative key rate, traced to detector clicks that did not appear
statistically coupled to the underlying phase/SNS-category variables in that
particular dataset — a property of that dataset, not a defect in the
pipeline (the report's own recommended next step is validating the pipeline
against a properly noise-calibrated simulation run).

**Expected contents** (add the actual pipeline files here):

```
preprocessing_pipeline/
├── pipeline.ipynb / pipeline.py   # the Task 2 notebook or script
├── inputs/                         # the 5 CSVs consumed (from Task 1)
│   ├── alice_signal_record.csv
│   ├── bob_signal_record.csv
│   ├── charlie_record_x.csv
│   ├── charlie_record_z.csv
│   └── charlie_xbasis_log.csv
└── outputs/                        # computed tables / figures
```

**The pipeline notebook/script and its CSV inputs were not included when
this folder structure was set up — add them here directly.**

## Why this is kept separate from `../qsafire`

Both `osd/` and `qsafire/` compute a "secure key rate," but from unrelated
channel models: a 40 km fibre link here, versus several-hundred-to-
1,000+ km satellite-to-ground free-space links in `qsafire/`. No file, CSV
schema, or code path is shared between the two — see
`../Summer_Internship_Report_DRDO.pdf`, Section 1.2 and Figure 1.1, for the
report's own statement of this boundary.
