"""
Water-mass boxes of the WHIRLS region in S_A-Theta space (TEOS-10), and the classification of
an MVP section. Definitions from
writing/Agulhas_Retroflection_water_masses_WHIRLS_FINAL_origins.docx (see anom_water_masses.ipynb
for the explanation of every limit).

  WM_BOXES                 list of boxes: id, name, group, cite, pub, lims, p (depth range), color
  WM_COLORS                colour of each box (surface yellows, central blues, mode greens,
                           intermediate pinks, deep greys)
  WM_BY_ID, WM_IDS         lookup by id, ids in order
  box_mask(b, SA, CT)      True where (S_A, Theta) is inside box b
  depth_mask(b, P)         True where pressure P is inside the depth range of box b
  classify(sec)            masks of every box on a section grid (S_A-Theta box AND depth range)
  box_text(b), depth_text(b)
"""

import gsw
import numpy as np


def _wm(id, name, group, cite, pub, theta=None, SP=None, CT=None, SA=None, sig0=None, sig1=None):
    """One box. Limits are (lo, hi), inclusive, None = open.
    theta : published potential temperature -> the curve Theta = gsw.CT_from_pt(S_A, theta)
    SP    : published practical salinity    -> S_A limits with gsw.SR_from_SP
    CT, SA, sig0 (= sigma_theta), sig1 : used as given (functions of S_A, Theta only)."""
    lims = {}
    if SP is not None:
        lims["SA"] = tuple(None if s is None else float(gsw.SR_from_SP(s)) for s in SP)
    for key, val in [("SA", SA), ("theta", theta), ("CT", CT), ("sig0", sig0), ("sig1", sig1)]:
        if val is not None:
            lims[key] = val
    return dict(id=id, name=name, group=group, cite=cite, pub=pub, lims=lims)

wm_transects = {"Transect 1": ["IOSTMW2", "AAIW_X22", "RSW_B00", "IEW", "ICW_SV42", "NICW_YT93", "AAMW_YT93"]}
EM86, BAR25 = "Emery and Meincke (1986)", "Barabinot et al. (2025)"
WM_BOXES = [
    # ---- surface waters ----
    _wm("SMW", "Salinity Maximum Water / Subtropical Underwater", "surface", BAR25, "σ0 < 25.6", sig0=(None, 25.6)),
    # _wm("SASW", "Subantarctic Surface Water", "surface", EM86, "θ 3.2–15.0, S 34.0–35.5", theta=(3.2, 15.0), SP=(34.0, 35.5)),
    _wm("AASW", "Antarctic Surface Water", "surface", EM86, "θ −1.0–1.0, S 34.0–34.6", theta=(-1.0, 1.0), SP=(34.0, 34.6)),
    _wm("BBW", "Bengal Bay Water", "surface", EM86, "θ 25–29, S 28.0–35.0", theta=(25.0, 29.0), SP=(28.0, 35.0)),
    _wm("ASW", "Arabian Sea Water", "surface", EM86, "θ 24–30, S 35.5–36.8", theta=(24.0, 30.0), SP=(35.5, 36.8)),
    _wm("IEW", "Indian Equatorial Water", "surface", EM86, "θ 8–23, S 34.6–35.0", theta=(8.0, 23.0), SP=(34.6, 35.0)),
    # _wm("IUW", "Indonesian Upper Water", "surface", EM86, "θ 8–23, S 34.4–35.0", theta=(8.0, 23.0), SP=(34.4, 35.0)),
    # ---- central (thermocline) waters ----
    _wm("SICW_EM86", "South Indian Central Water", "central", EM86, "θ 8–25, S 34.6–35.8", theta=(8.0, 25.0), SP=(34.6, 35.8)),
    _wm("ICW_SV42", "Indian Central Water", "central", "You (1997), after Sverdrup et al. (1942)",
        "θ 7–16, S 34.70–35.70", theta=(7.0, 16.0), SP=(34.70, 35.70)),
    _wm("NICW_YT93", "North Indian Central Water", "central", "You (1997), after You and Tomczak (1993)",
        "θ 7.80–15.72, S 34.84–35.10", theta=(7.80, 15.72), SP=(34.84, 35.10)),
    _wm("AAMW_YT93", "Australasian Mediterranean Water", "central", "You (1997), after You and Tomczak (1993)",
        "θ 5.55–14.89, S 34.52–34.77", theta=(5.55, 14.89), SP=(34.52, 34.77)),
    _wm("SACW", "South Atlantic Central Water", "central", BAR25, "σ0 25.6–26.5", sig0=(25.6, 26.5)),
    _wm("ESACW", "Eastern South Atlantic Central Water", "central", BAR25, "Θ 9.44–13.60, S_A 34.9–35.40, σ0 26.5–26.93",
        CT=(9.44, 13.60), SA=(34.90, 35.40), sig0=(26.5, 26.93)),
    # ---- mode waters ----
    _wm("SIMW", "South Indian Subtropical Mode Water", "mode", "Olson et al. (1992)", "θ 17.4–17.8", theta=(17.4, 17.8)),
    _wm("IOSTMW1", "Indian Ocean Subtropical Mode Water", "mode", "Tsubouchi et al. (2010)",
        "θ 16.54±0.49, S 35.51±0.04, σθ 26.0±0.1", theta=(16.05, 17.03), SP=(35.47, 35.55), sig0=(25.9, 26.1)),
    _wm("IOSTMW2", "Indian Ocean Subtropical Mode Water type 2", "mode", "Tsubouchi et al. (2010)",
        "θ 12.84±0.57, S 35.17±0.11, σθ 26.57±0.04", theta=(12.27, 13.41), SP=(35.06, 35.28), sig0=(26.53, 26.61)),
    _wm("SAMW_KT97", "Subantarctic Mode Water, SE Indian Ocean", "mode", "Karstensen and Tomczak (1997)",
        "θ 8.5–9.5, σθ 26.8–26.9", theta=(8.5, 9.5), sig0=(26.8, 26.9)),
    # ---- intermediate waters ----
    _wm("AAIW_IO", "Antarctic Intermediate Water, Indian Ocean", "intermediate", EM86, "θ 2–10, S 33.8–34.8",
        theta=(2.0, 10.0), SP=(33.8, 34.8)),
    # _wm("AAIW_AO", "Antarctic Intermediate Water, Atlantic", "intermediate", EM86, "θ 2–6, S 33.8–34.8",
    #     theta=(2.0, 6.0), SP=(33.8, 34.8)),
    # _wm("AAIW_LT21", "Antarctic Intermediate Water", "intermediate", BAR25 + ", after Liu and Tanhua (2021)",
    #     "σ0 26.9–32.15", sig0=(26.9, 32.15)),
    _wm("AAIW_X22", "AAIW, salinity-minimum layer", "intermediate", BAR25 + ", after Xia et al. (2022)",
        "σ0 27.0–27.2", sig0=(27.0, 27.2)),
    _wm("IIW", "Indonesian Intermediate Water", "intermediate", EM86, "θ 3.5–5.5, S 34.6–34.7",
        theta=(3.5, 5.5), SP=(34.6, 34.7)),
    _wm("RSPGIW", "Red Sea–Persian Gulf Intermediate Water", "intermediate", EM86, "θ 5–14, S 34.8–35.4",
        theta=(5.0, 14.0), SP=(34.8, 35.4)),
    _wm("RSW_B00", "Red Sea Water layer", "intermediate", "Beal et al. (2000)", "σ1 31.5–32.1", sig1=(31.5, 32.1)),
    # ---- deep and abyssal waters (below MVP depth, CTD only) ----
    # _wm("NADW", "North Atlantic Deep Water", "deep", EM86, "θ 1.5–4.0, S 34.8–35.0", theta=(1.5, 4.0), SP=(34.8, 35.0)),
    # _wm("CDW", "Circumpolar Deep Water", "deep", EM86, "θ 0.1–2.0, S 34.62–34.73", theta=(0.1, 2.0), SP=(34.62, 34.73)),
    # _wm("AABW", "Antarctic Bottom Water", "deep", EM86, "θ −0.9–1.7, S 34.64–34.72", theta=(-0.9, 1.7), SP=(34.64, 34.72)),
]
# ---- depth check: a box only counts inside the pressure range of its group ----
# A sample in the S_A-Theta box but outside this range is REMOVED from the box (and reported).
# SACW and ESACW are exempt: the Agulhas region is where they form, so they are not yet tied to a depth.
DEPTH_RANGES = {"surface": (0, 250), "central": (100, 800), "mode": (0, 800),
                "intermediate": (250, 2000), "deep": (1500, None)}        # dbar, inclusive, None = open
DEPTH_EXEMPT = {"SACW", "ESACW"}
for b in WM_BOXES:
    b["p"] = None if b["id"] in DEPTH_EXEMPT else DEPTH_RANGES[b["group"]]

# ---- one colour per water mass, shades of one hue per group (light -> dark within the group):
# surface = yellows, central = blues, mode = greens, intermediate = pinks, deep = greys.
# Ids that are commented out above have a colour too, so they can be switched back on.
WM_COLORS = {
    # surface: yellows
    "SMW": "#FFF3A6", "SASW": "#FFE66B", "AASW": "#F5D33A", "BBW": "#E3BC12",
    "ASW": "#C9A300", "IEW": "#A68700", "IUW": "#7F6800",
    # central: blues
    "SICW_EM86": "#C6DBEF", "ICW_SV42": "#9ECAE1", "NICW_YT93": "#6BAED6",
    "AAMW_YT93": "#4292C6", "SACW": "#2171B5", "ESACW": "#08519C",
    # mode: greens
    "SIMW": "#C7E9C0", "IOSTMW1": "#74C476", "IOSTMW2": "#31A354", "SAMW_KT97": "#006D2C",
    # intermediate: pinks
    "AAIW_IO": "#FDE0EF", "AAIW_AO": "#F9C4E1", "AAIW_LT21": "#F49AC9", "AAIW_X22": "#E46FAE",
    "IIW": "#D13F8F", "RSPGIW": "#B01870", "RSW_B00": "#7A0150",
    # deep: greys
    "NADW": "#BDBDBD", "CDW": "#737373", "AABW": "#252525",
}
GROUP_DEFAULT_COLORS = {"surface": "#F5D33A", "central": "#4292C6", "mode": "#31A354",
                        "intermediate": "#D13F8F", "deep": "#737373"}
for b in WM_BOXES:                     # a new box without an entry gets its group's colour
    b["color"] = WM_COLORS.get(b["id"], GROUP_DEFAULT_COLORS[b["group"]])

WM_IDS = [b["id"] for b in WM_BOXES]
WM_BY_ID = {b["id"]: b for b in WM_BOXES}
assert len(set(WM_IDS)) == len(WM_IDS) <= 62        # one bit per box in an int64 code


def box_mask(b, SA, CT, sig0=None, sig1=None):
    """True where (S_A, Theta) satisfies every limit of box b (limits inclusive)."""
    SA, CT = np.asarray(SA, float), np.asarray(CT, float)
    m = np.isfinite(SA) & np.isfinite(CT)
    with np.errstate(invalid="ignore"):
        for var, (lo, hi) in b["lims"].items():
            if var == "theta":                    # theta limit = curve Theta(S_A) in the S_A-Theta plane
                if lo is not None:
                    m &= CT >= gsw.CT_from_pt(SA, lo)
                if hi is not None:
                    m &= CT <= gsw.CT_from_pt(SA, hi)
                continue
            if var == "sig0":
                x = gsw.sigma0(SA, CT) if sig0 is None else sig0
            elif var == "sig1":
                x = gsw.sigma1(SA, CT) if sig1 is None else sig1
            else:
                x = {"SA": SA, "CT": CT}[var]
            if lo is not None:
                m &= x >= lo
            if hi is not None:
                m &= x <= hi
    return m


def depth_mask(b, P):
    """True where pressure P (dbar) is inside the depth range of box b (always True if exempt)."""
    P = np.asarray(P, float)
    if b["p"] is None:
        return np.ones(P.shape, bool)
    lo, hi = b["p"]
    with np.errstate(invalid="ignore"):
        return (P >= (lo if lo is not None else -np.inf)) & (P <= (hi if hi is not None else np.inf))


def depth_text(b):
    if b["p"] is None:
        return "no depth limit"
    lo, hi = b["p"]
    return f"p ≥ {lo:g} dbar" if hi is None else f"p {lo:g}–{hi:g} dbar"


def box_text(b):
    """Box in its S_A-Theta form, e.g. 'Θ≈12.30–13.45 °C, S_A 35.23–35.45, σ0 26.53–26.61'.
    A theta limit is a curve; Θ≈ gives it at the middle of the box's S_A range (35 if none)."""
    L = b["lims"]

    def rng(lo, hi, f="{:g}"):
        if lo is None:
            return "≤ " + f.format(hi)
        if hi is None:
            return "≥ " + f.format(lo)
        return f.format(lo) + "–" + f.format(hi)

    sa = [s for s in L.get("SA", (None, None)) if s is not None]
    sa_mid = float(np.mean(sa)) if sa else 35.0
    parts = []
    if "theta" in L:
        lo, hi = (None if t is None else float(gsw.CT_from_pt(sa_mid, t)) for t in L["theta"])
        parts.append("Θ≈" + rng(lo, hi, "{:.2f}") + " °C")
    if "CT" in L:
        parts.append("Θ " + rng(*L["CT"]) + " °C")
    if "SA" in L:
        parts.append("S_A " + rng(*L["SA"], "{:.2f}"))
    for k, s in [("sig0", "σ0"), ("sig1", "σ1")]:
        if k in L:
            parts.append(f"{s} " + rng(*L[k]))
    return ", ".join(parts)


def classify(sec):
    """Masks of every box on the section grid (S_A-Theta box AND depth range of its group),
    the samples removed by the depth check, and an int64 code per sample (bit k = in box k)."""
    SA, CT = sec["SA"], sec["CT"]
    P = np.broadcast_to(sec["p"][:, None], SA.shape)
    sig0, sig1 = gsw.sigma0(SA, CT), gsw.sigma1(SA, CT)
    masks, removed = {}, {}
    for b in WM_BOXES:
        in_ts = box_mask(b, SA, CT, sig0, sig1)
        in_p = depth_mask(b, P)
        masks[b["id"]] = in_ts & in_p
        removed[b["id"]] = in_ts & ~in_p
    code = np.zeros(SA.shape, np.int64)
    for k, i in enumerate(WM_IDS):
        code |= masks[i].astype(np.int64) << k
    return masks, removed, code, np.isfinite(SA) & np.isfinite(CT)


def combo_ids(code):
    return [i for k, i in enumerate(WM_IDS) if (int(code) >> k) & 1]


def combo_name(code, sep=" + "):
    return sep.join(combo_ids(code)) or "no box"
