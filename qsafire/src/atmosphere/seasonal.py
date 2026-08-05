"""Seasonal aerosol/turbulence dataset (notebook Part VI: explicit winter /
summer / monsoon seasonal extension of Phases 3-4).

Every value is traceable to a specific cited source (see notebook Section
6.1 / cells 132-134, 136 for full citations: R2007, L2013, SB2011, G2006,
A2013, GD2018). Confidence tags follow the High/Medium/Indicative convention
established in Phase 3/4.

Origin: notebook cells 135 (SEASONAL_AEROSOL), 137 (A0_BASELINE,
SEASONAL_TURBULENCE).
"""

# Seasonal AOD(500nm)/Angstrom parameter sets -- see notebook Section 6.3 for citations
SEASONAL_AEROSOL = {
    "Delhi (Alice)": {
        "Winter (DJF)":              dict(aod0=0.70, alpha=0.90, lam0_nm=500.0, conf="Medium (spread across R2007/G2006/SB2011, see 6.3)"),
        "Summer/pre-monsoon (MAMJ)": dict(aod0=1.00, alpha=0.50, lam0_nm=500.0, conf="High (R2007 June peak; SB2011 seasonal alpha minimum)"),
        "Monsoon (JAS)":             dict(aod0=0.90, alpha=0.70, lam0_nm=500.0, conf="Medium (regional monsoon AOD; L2013 qualitative coarse-mode confirmation)"),
    },
    "Mumbai (Bob)": {
        "Winter (DJF)":              dict(aod0=0.31, alpha=0.90, lam0_nm=550.0, conf="High (R2007 Results, direct multi-year mean)"),
        "Summer/pre-monsoon (MAMJ)": dict(aod0=0.74, alpha=0.90, lam0_nm=550.0, conf="High (R2007 Results, July peak)"),
        "Monsoon (JAS)":             dict(aod0=0.65, alpha=0.90, lam0_nm=550.0, conf="Indicative (interpolated, no JAS-only mean in R2007)"),
    },
}

# Seasonal ground-level turbulence scaling factors -- see notebook Section 6.4 for citations
A0_BASELINE = 1.7e-14
SEASONAL_TURBULENCE = {
    "Delhi (Alice)": {
        "Winter (DJF)":              dict(A0=A0_BASELINE*7.08, conf="Medium-proxy (A2013 radio-Cn2 seasonal ratio, see 6.2/6.4)"),
        "Summer/pre-monsoon (MAMJ)": dict(A0=A0_BASELINE*0.63, conf="Medium-proxy (A2013 radio-Cn2 seasonal ratio)"),
        "Monsoon (JAS)":             dict(A0=A0_BASELINE*1.00, conf="Medium-proxy (A2013 reference season)"),
    },
    "Mumbai (Bob)": {
        "Winter (DJF)":              dict(A0=A0_BASELINE*1.00, conf="Indicative (GD2018 qualitative only)"),
        "Summer/pre-monsoon (MAMJ)": dict(A0=A0_BASELINE*1.00, conf="Indicative (GD2018 qualitative only)"),
        "Monsoon (JAS)":             dict(A0=A0_BASELINE*0.50, conf="Indicative (GD2018: monsoon near-surface minimum)"),
    },
}

SEASONS = ["Winter (DJF)", "Summer/pre-monsoon (MAMJ)", "Monsoon (JAS)"]
