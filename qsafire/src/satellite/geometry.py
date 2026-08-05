"""Ground-to-satellite geometry engine (notebook Part III / Phase 2).

Standard textbook astrodynamics/geodesy, no protocol-specific assumptions:
WGS-84 geodetic <-> ECEF, sun-synchronous circular-orbit design (solved, not
assumed, via the J2 secular nodal-precession condition), two-body + first-order
J2 secular propagation, GMST Earth orientation, and ENU topocentric look
angles (elevation/azimuth/slant-range).

Origin: notebook cells 80 (constants + SITES), 82 (geodetic_to_ecef,
enu_rotation), 84 (sun_sync_orbit), 87 (julian_date, gmst_rad, eci_position,
eci_to_ecef, look_angles, refraction_bump_deg), 89 (slant_from_elev), 91
(gc_midpoint, sub_satellite, haversine_km).
"""

import numpy as np

# ============================================================
# 3.0  Constants (WGS-84 geodesy + standard gravitational parameters)
# ============================================================
MU_EARTH = 398600.4418            # km^3/s^2  (GM, geodetic standard, e.g. Vallado App. D)
RE_EQ    = 6378.137                # km, WGS-84 equatorial radius
F_WGS84  = 1/298.257223563         # WGS-84 flattening
E2_WGS84 = 2*F_WGS84 - F_WGS84**2  # first eccentricity squared
J2       = 1.08262668e-3           # Earth's J2 zonal harmonic (dimensionless)
OMEGA_E  = 7.2921150e-5            # rad/s, Earth sidereal rotation rate
YEAR_TROPICAL_DAYS = 365.2421897   # days, mean tropical year

# Ground stations (Phase 1 document, Section 3.1) -- WGS-84 city reference points
SITES = {
    "Delhi (Alice)": dict(lat=28.6139, lon=77.2090, alt_km=0.216),
    "Mumbai (Bob)":  dict(lat=19.0760, lon=72.8777, alt_km=0.014),
}


def geodetic_to_ecef(lat_deg, lon_deg, alt_km):
    """WGS-84 geodetic (lat, lon, altitude) -> ECEF cartesian, km. Cell 82."""
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    N = RE_EQ / np.sqrt(1 - E2_WGS84*np.sin(lat)**2)
    x = (N + alt_km)*np.cos(lat)*np.cos(lon)
    y = (N + alt_km)*np.cos(lat)*np.sin(lon)
    z = (N*(1-E2_WGS84) + alt_km)*np.sin(lat)
    return np.array([x, y, z])


def enu_rotation(lat_deg, lon_deg):
    """Rows are the East, North, Up unit vectors expressed in ECEF, at (lat, lon). Cell 82."""
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    east  = np.array([-np.sin(lon),             np.cos(lon),            0.0])
    north = np.array([-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)])
    up    = np.array([ np.cos(lat)*np.cos(lon),  np.cos(lat)*np.sin(lon), np.sin(lat)])
    return np.vstack([east, north, up])


def build_site_ecef():
    """ECEF position vectors for every entry in SITES. Cell 82 (script part)."""
    return {name: geodetic_to_ecef(s['lat'], s['lon'], s['alt_km']) for name, s in SITES.items()}


def sun_sync_orbit(h_km, Re=RE_EQ, mu=MU_EARTH, j2=J2):
    """Circular sun-synchronous orbit: solves for inclination, returns
    first-order J2 secular rates. Cell 84."""
    a = Re + h_km
    n0 = np.sqrt(mu/a**3)
    p = a  # circular orbit, e=0
    omega_sun = 2*np.pi/(YEAR_TROPICAL_DAYS*86400.0)
    cos_i = -omega_sun/(1.5*n0*j2*(Re/p)**2)
    if abs(cos_i) > 1:
        raise ValueError(f"No sun-synchronous circular solution at h={h_km} km with this model")
    i = np.arccos(cos_i)
    raan_dot = -1.5*n0*j2*(Re/p)**2*np.cos(i)
    argp_dot =  1.5*n0*j2*(Re/p)**2*(2 - 2.5*np.sin(i)**2)
    n_true   = n0*(1 + 1.5*j2*(Re/p)**2*(1 - 1.5*np.sin(i)**2))
    u_dot    = n_true + argp_dot
    return dict(a=a, i=i, n0=n0, n_true=n_true, raan_dot=raan_dot, u_dot=u_dot)


def julian_date(year, month, day, hour=0, minute=0, sec=0.0):
    """Cell 87."""
    if month <= 2:
        year -= 1; month += 12
    A = year // 100
    B = 2 - A + A//4
    jd0 = int(365.25*(year+4716)) + int(30.6001*(month+1)) + day + B - 1524.5
    return jd0 + (hour + minute/60 + sec/3600)/24.0


def gmst_rad(jd):
    """IAU-1982-style low-precision GMST (deg->rad), standard textbook series. Cell 87."""
    T = (jd - 2451545.0)/36525.0
    gmst_deg = (280.46061837 + 360.98564736629*(jd - 2451545.0)
                + 0.000387933*T**2 - T**3/38710000.0)
    return np.radians(gmst_deg % 360.0)


def eci_position(a, i, raan, u):
    """Cell 87."""
    r_orb = a*np.array([np.cos(u), np.sin(u), 0.0])
    Rz = np.array([[np.cos(raan), -np.sin(raan), 0], [np.sin(raan), np.cos(raan), 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, np.cos(i), -np.sin(i)], [0, np.sin(i), np.cos(i)]])
    return Rz @ Rx @ r_orb


def eci_to_ecef(r_eci, gmst):
    """Cell 87."""
    Rz = np.array([[np.cos(gmst), np.sin(gmst), 0], [-np.sin(gmst), np.cos(gmst), 0], [0, 0, 1]])
    return Rz @ r_eci


def look_angles(sat_ecef, site_ecef_vec, site_lat, site_lon):
    """Cell 87."""
    rho = sat_ecef - site_ecef_vec
    R = enu_rotation(site_lat, site_lon)
    e, n, u = R @ rho
    rng = np.linalg.norm(rho)
    el = np.degrees(np.arcsin(u/rng))
    az = np.degrees(np.arctan2(e, n)) % 360
    return el, az, rng


def refraction_bump_deg(true_elev_deg):
    """Bennett (1982) approximate atmospheric refraction (arcmin), for
    reference/diagnostics. Cell 87."""
    if true_elev_deg < -1:
        return 0.0
    R_arcmin = 1.0/np.tan(np.radians(true_elev_deg + 7.31/(true_elev_deg + 4.4)))
    return R_arcmin/60.0


def slant_from_elev(elev_deg, h_km, Re=RE_EQ):
    """Spherical-Earth elevation <-> slant-range relation. Cell 89."""
    eps = np.radians(elev_deg)
    ratio = (Re + h_km)/Re
    return Re*(np.sqrt(ratio**2 - np.cos(eps)**2) - np.sin(eps))


def gc_midpoint(lat1, lon1, lat2, lon2):
    """Great-circle midpoint of two (lat,lon) points, degrees in/out. Cell 91."""
    p1, p2 = np.radians([lat1, lon1]), np.radians([lat2, lon2])
    dlon = p2[1]-p1[1]
    Bx, By = np.cos(p2[0])*np.cos(dlon), np.cos(p2[0])*np.sin(dlon)
    lat_m = np.arctan2(np.sin(p1[0])+np.sin(p2[0]), np.sqrt((np.cos(p1[0])+Bx)**2 + By**2))
    lon_m = p1[1] + np.arctan2(By, np.cos(p1[0])+Bx)
    return np.degrees(lat_m), np.degrees(lon_m)


def sub_satellite(t, a, inc, raan0, u0, raan_dot, u_dot, gmst0):
    """Cell 91."""
    raan = raan0 + raan_dot*t
    u = u0 + u_dot*t
    lat = np.arcsin(np.sin(inc)*np.sin(u))
    lon = raan + np.arctan2(np.cos(inc)*np.sin(u), np.cos(u)) - (gmst0 + OMEGA_E*t)
    return np.degrees(lat), (np.degrees(lon)+180) % 360 - 180


def haversine_km(lat1, lon1, lat2, lon2, R=RE_EQ):
    """Cell 91."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2-lat1); dl = np.radians(lon2-lon1)
    x = np.sin(dphi/2)**2 + np.cos(p1)*np.cos(p2)*np.sin(dl/2)**2
    return 2*R*np.arcsin(np.sqrt(x))
