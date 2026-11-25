"""
Water-mass definitions collected from the literature (literature review, WHIRLS).

One entry per (water mass, study). T is potential temperature theta (°C) and S is practical
salinity (psu), as in the source tables, so the comparison is done in (theta, S_P) and not in
(Theta, S_A).

Fields
  wm       short name of the water mass (one colour per wm in the figures)
  name     full name
  study    the reference that gives the definition
  T, S     (min, max) range; None for an open bound ("> 7" -> (7, None))
  T0, S0   single core value, when the study gives a point rather than a range
  sig      optional (min, max) sigma0 range
  notes    free text (how the study defines it)
  table    where the definition was copied from (the table that cites the study)

Entries with a range (T and S) are used to flag the section; entries with only a core value
(T0, S0) appear only as stars on the theta-S diagram; entries with neither are kept for the
record (notes only).
"""

TABLE_INDIAN = "Indian Ocean water-mass table (screenshot sent 2026-11-17)"

WM_NAMES = {
    "ICW": "Indian Central Water",
    "AAMW": "Australasian Mediterranean Water",
    "NICW": "North Indian Central Water",
    "RSW": "Red Sea Water",
    "BBW": "Bengal Bay Water",
    "ASW": "Arabian Sea Water",
    "IEW": "Indian Equatorial Water",
}


def _e(wm, study, T=None, S=None, T0=None, S0=None, sig=None, notes="", table=TABLE_INDIAN):
    return dict(wm=wm, name=WM_NAMES[wm], study=study, T=T, S=S, T0=T0, S0=S0, sig=sig,
                notes=notes, table=table)


WM_LIT = [
    # Indian Central Water
    _e("ICW", "Sverdrup et al. (1942)", T=(7.00, 16.00), S=(34.70, 35.70)),
    _e("ICW", "Mamayev (1975)", T0=16.00, S0=35.60),
    _e("ICW", "Emery and Meincke (1986)", T=(8.00, 25.00), S=(34.60, 35.80)),
    _e("ICW", "Warren et al. (1966)", T=(7.00, 15.00), S=(34.70, 35.50)),
    _e("ICW", "McCartney (1977)", notes="Formed by late winter convective overturning at 40-45°S."),
    # Australasian Mediterranean Water
    _e("AAMW", "You and Tomczak (1993)", T=(5.55, 14.89), S=(34.52, 34.77)),
    _e("AAMW", "Mamayev (1975)", T0=25.00, S0=34.50),
    _e("AAMW", "Emery and Meincke (1986)", T=(8.00, 23.00), S=(34.40, 35.00)),
    _e("AAMW", "Wyrtki (1961)", notes="After Warren (1981)."),
    _e("AAMW", "Rochford (1966)"),
    # North Indian Central Water
    _e("NICW", "Gordon (1986)", notes="Thermocline water below the surface water in the Bay of Bengal."),
    _e("NICW", "You and Tomczak (1993)", T=(7.80, 15.72), S=(34.84, 35.10),
       notes="Same as Gordon (1986) but defined as aged ICW."),
    # Red Sea Water
    _e("RSW", "Sverdrup et al. (1942)", T=(7.00, None), S=(35.00, None),
       notes="Defined as an oxygen minimum core (T/S bounds alone are not specific)."),
    _e("RSW", "Mamayev (1975)", T0=23.00, S0=40.00, notes="Defined as intermediate water."),
    # Bengal Bay Water
    _e("BBW", "Mamayev (1975)", T0=25.00, S0=33.80, notes="Surface low-salinity water."),
    _e("BBW", "Emery and Meincke (1986)", T=(25.00, 29.00), S=(28.00, 35.00),
       notes="Surface low-salinity water."),
    # Arabian Sea Water
    _e("ASW", "Emery and Meincke (1986)", T=(24.00, 30.00), S=(35.50, 36.80)),
    # Indian Equatorial Water
    _e("IEW", "Sverdrup et al. (1942)", T=(5.00, 16.00), S=(34.90, 35.20)),
    _e("IEW", "Mamayev (1975)", T0=25.00, S0=35.30),
    _e("IEW", "Emery and Meincke (1986)", T=(8.00, 23.00), S=(34.60, 35.00)),
]

# Okabe-Ito (colour-blind safe), one colour per water mass, in the order of WM_NAMES
WM_LIT_COLORS = dict(zip(WM_NAMES, ["#E69F00", "#56B4E9", "#009E73", "#D55E00",
                                    "#0072B2", "#CC79A7", "#000000"]))


def in_range(val, bounds):
    lo, hi = bounds
    out = val == val                      # False where NaN
    if lo is not None:
        out &= val >= lo
    if hi is not None:
        out &= val <= hi
    return out


def lit_mask(entry, pt, SP, sig0=None):
    """True where (theta, S_P[, sigma0]) falls inside the entry's ranges."""
    m = in_range(pt, entry["T"]) & in_range(SP, entry["S"])
    if entry["sig"] is not None and sig0 is not None:
        m &= in_range(sig0, entry["sig"])
    return m


def ranged(wm=None):
    """Entries that have a T and S range (optionally for one water mass)."""
    return [e for e in WM_LIT if e["T"] is not None and e["S"] is not None
            and (wm is None or e["wm"] == wm)]
