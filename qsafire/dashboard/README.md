# QSAFire Dashboard

An interactive Streamlit dashboard over the full satellite-link pipeline
(Parts III-XIII). Run:

```bash
streamlit run dashboard/app.py
```

then open the URL Streamlit prints (typically `http://localhost:8501`).

## What it shows

| Tab | Content |
|---|---|
| Pass Geometry | Delhi-Mumbai elevation/slant-range vs. time, the pass summary table, ground-site map |
| Atmosphere | Zenith optical depth table, atmospheric transmittance vs. time |
| Turbulence | Rytov variance vs. time (shared baseline vs. site case), Gamma-Gamma regime mix |
| Real-Loss Key Rate | Fixed-decoy vs. joint-optimized key rate vs. time, root-cause loss breakdown, key rate vs. elevation |
| Weather & Annual | Weather-gated annual expected key volume, old-vs-corrected seasonal robustness, sensitivity table |
| Multi-City Links | The four evaluated links' summary table, bottleneck key-rate comparison, 6-city map |
| Final Metrics | The fully asymmetric per-link master table (QBER, phase error, yield, SKR) and comparison bar charts |

Headline metrics (Delhi-Mumbai distance/SKR, best link, common-visibility
window duration) are pinned at the top of every tab.

## Performance

The full pipeline runs Monte Carlo ensemble optimization across many points
in the pass and four city-pair links, so the **first load takes on the order
of a few minutes**. The result is cached for the rest of the session
(`st.cache_data`) — reopening or switching tabs afterward is instant. Use the
"Run / refresh full simulation" button in the sidebar to force a recompute
(e.g. after editing a module's constants).

## Requirements

`streamlit` (added to `requirements.txt` / `pyproject.toml`). Everything
else the dashboard imports is already part of the package
(`src.experiments.*`, `src.satellite.geometry`, `src.links.city_pairs`).
