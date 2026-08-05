"""City parameter database for the four evaluated links (notebook Part XII /
Phase 9). Delhi/Mumbai repeated here for a uniform single-table comparison;
all values/sources as established in Parts III-VI or newly cited in notebook
Section 12.1 (Sagar et al./ARIES for Nainital, Ramachandran & Kedia 2011/2012
for Ahmedabad, Fenner et al. ACP 2022 + Atmos. Environ. 2013 for Bangalore,
Kaskaoutis et al. 2009 + a 2020s MODIS trend study for Hyderabad).

Origin: notebook cell 227.
"""

CITY_DB = {
    "Delhi":      dict(lat=28.6139, lon=77.2090, alt_km=0.216, aod500=0.78, alpha=0.78,
                        turb_mult=1.0, conf="High (Soni et al. 2011)"),
    "Mumbai":     dict(lat=19.0760, lon=72.8777, alt_km=0.014, aod500=0.45, alpha=0.90,
                        turb_mult=1.0, conf="High (Ramachandran 2007)"),
    "Nainital":   dict(lat=29.3667, lon=79.4586, alt_km=1.951, aod500=0.18, alpha=1.00,
                        turb_mult=0.4, conf="Medium (Sagar et al., ARIES report; site-selection proxy for turbulence)"),
    "Ahmedabad":  dict(lat=23.0225, lon=72.5714, alt_km=0.053, aod500=0.40, alpha=1.00,
                        turb_mult=1.0, conf="High (Ramachandran & Kedia 2011/2012)"),
    "Bangalore":  dict(lat=12.9716, lon=77.5946, alt_km=0.920, aod500=0.25, alpha=1.00,
                        turb_mult=1.0, conf="High (ACP 2022 megacity comparison; Atmos. Environ. 2013)"),
    "Hyderabad":  dict(lat=17.3850, lon=78.4867, alt_km=0.542, aod500=0.55, alpha=0.90,
                        turb_mult=1.0, conf="Medium (Kaskaoutis et al. 2009, bridged with a 2020s trend study)"),
}

LINK_PAIRS = [
    ("Delhi", "Mumbai"),
    ("Delhi", "Nainital"),
    ("Mumbai", "Ahmedabad"),
    ("Bangalore", "Hyderabad"),
]

LINK_COLORS = {"Delhi-Mumbai": "#1f77b4", "Delhi-Nainital": "#2ca02c",
               "Mumbai-Ahmedabad": "#d62728", "Bangalore-Hyderabad": "#9467bd"}
