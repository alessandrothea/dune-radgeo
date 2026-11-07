"""Detector profiles: volume name patterns for VD and HD geometries."""
from __future__ import annotations

_VD_STRUCT = [r"^volTPCActive$", r"^volGaseousArgon$", r"^volTPC$", r"^volCryostat$"]

PROFILES: dict[str, dict] = {
    # Profiles are per detector and GDML flavour (perl- or ggd-generated);
    # "detector" is what goes into the IR.  VD, ggd GDML: same detector volume
    # names as the perl GDML, but the cryostat sits in
    # volSupportEnc/volEnclosureCryostat, surrounded by passive material.
    # Checked before "vd-perl", whose hint volume it also contains.
    "vd-ggd": {
        "label": "VD (Vertical Drift), ggd GDML",
        "detector": "vd",
        "anode_patterns":  [r"^volAnodePlate$", r"^volAnodePlateBottom$"],
        "opdet_patterns":  [r"^volArapuca"],
        "cathode_patterns": [r"^volCathodeGrid$"],
        "fieldshaper_patterns": [r"^volFieldShaper(Slim)?$"],
        "struct_patterns": _VD_STRUCT,
        # --with-passive: cryostat insulation, support structure, shielding, cavern
        "passive_patterns": [r"^volSteelShell$", r"^fShellLog$", r"^WoodLog$",
                             r"^FoamLog$", r"^ShellOutLog$",
                             r"^IBeam(Top|Bot|Side)$", r"^Belt(Uni|UniTop|HoleUni)$",
                             r"^ShieldBlockContainerLog$", r"^boxshape(Lead)?Vol$",
                             r"^volShotbox$", r"^volGrout$", r"^volConcrete$",
                             r"^volRadioRockShell$"],
        "hint_vols": {"volSupportEnc", "volEnclosureCryostat"},
    },
    # VD, perl GDML
    "vd-perl": {
        "label": "VD (Vertical Drift), perl GDML",
        "detector": "vd",
        "anode_patterns":  [r"^volAnodePlate$", r"^volAnodePlateBottom$"],
        "opdet_patterns":  [r"^volArapuca"],
        # v5+ name; v4 geometries call the same CathodeGrid solid volGroundGrid
        "cathode_patterns": [r"^volCathodeGrid$", r"^volGroundGrid$"],
        "fieldshaper_patterns": [r"^volFieldShaper(Slim)?$"],
        "struct_patterns": _VD_STRUCT,
        "passive_patterns": [r"^volSteelShell$", r"^volFoamPadding$",
                             r"^volSteelSupport$"],
        "hint_vols": {"volAnodePlate"},
    },
    "hd": {
        "label": "HD (Horizontal Drift / APA)",
        "anode_patterns":  [r"^volAPAFrameYSide$"],
        "opdet_patterns":  [r"^volArapuca_\d"],
        "cathode_patterns": [],   # no HD reference GDML checked yet
        "fieldshaper_patterns": [],
        "struct_patterns": [r"^volTPCActiveInner$", r"^volGaseousArgon$",
                            r"^volTPC$", r"^volCryostat$"],
        "hint_vols": {"volAPAFrameYSide"},
    },
}

_GENERIC_PROFILE: dict = {
    "label": "DUNE (generic)",
    "anode_patterns":  [r"volAnode|volCRM$|APAFrame"],
    "opdet_patterns":  [r"volOpDetSensitive|volArapuca"],
    "cathode_patterns": [r"^volCathodeGrid$", r"^volGroundGrid$"],
    "fieldshaper_patterns": [r"^volFieldShaper(Slim)?$"],
    "struct_patterns": [r"volTPCActive$", r"volGaseousArgon$",
                        r"^volTPC$", r"^volCryostat$"],
}


def detect_profile(vol_names: set[str]) -> str | None:
    """First profile (in PROFILES order) with a hint volume in `vol_names`."""
    for key, prof in PROFILES.items():
        if prof["hint_vols"] & vol_names:
            return key
    return None


def profile_detector(key: str | None) -> str:
    """Detector family written to the IR ("vd" / "hd" / ""), e.g. vd-ggd → vd."""
    if key and key in PROFILES:
        return PROFILES[key].get("detector", key)
    return ""


def get_profile(key: str | None) -> tuple[str, dict]:
    """Return (key, profile_dict); falls back to generic if key is None."""
    if key and key in PROFILES:
        return key, PROFILES[key]
    return "", _GENERIC_PROFILE
