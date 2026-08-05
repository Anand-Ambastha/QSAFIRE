"""Weather climatology dataset & availability proxy (notebook Part IX / Phase 7).

Real IMD station-climatology rainy-day counts (Safdarjung/Delhi, Santacruz/
Mumbai) plus a real Indo-Gangetic-Plain winter-fog climatology, combined into
an explicitly-flagged (Medium confidence) proxy for cloud-free-line-of-sight
probability -- see notebook Section 9.1 for full citations ([IMD-D], [IMD-M],
[MAUSAM-Fog]).

Origin: notebook cell 180.
"""

RAINY_DAYS_PER_MONTH = {
    "Delhi (Alice)": {1: 2, 2: 2, 3: 1, 4: 1, 5: 2, 6: 4, 7: 13, 8: 19, 9: 7, 10: 1, 11: 0, 12: 1},
    "Mumbai (Bob)":  {1: 0.1, 2: 0.1, 3: 0.1, 4: 0.0, 5: 0.6, 6: 14.1, 7: 22.1, 8: 20.2, 9: 14.0, 10: 3.6, 11: 0.5, 12: 0.3},
}
DAYS_IN_MONTH = {1: 31, 2: 28, 3: 31, 4: 30, 5: 31, 6: 30, 7: 31, 8: 31, 9: 30, 10: 31, 11: 30, 12: 31}

# Delhi-specific additional winter fog outage (Section 9.1, MAUSAM-Fog), applied Dec/Jan on top of
# (not double-counted with) rainy-day fraction -- fog days assumed distinct from rainy days.
FOG_DAYS_DELHI = {12: 8.5, 1: 8.0}  # representative of the 7-10 (Dec) / >8 (Jan) IGP range


def availability_fraction(station, month):
    """Cell 180."""
    rainy = RAINY_DAYS_PER_MONTH[station][month]
    fog = FOG_DAYS_DELHI.get(month, 0.0) if "Delhi" in station else 0.0
    outage_days = min(rainy + fog, DAYS_IN_MONTH[month])  # cap at 100% outage
    return 1.0 - outage_days/DAYS_IN_MONTH[month]
