"""Detector profiles: volume name patterns for VD and HD geometries."""
from __future__ import annotations

PROFILES: dict[str, dict] = {
    "vd": {
        "label": "VD (Vertical Drift)",
        "anode_patterns":  [r"^volAnodePlate$", r"^volAnodePlateBottom$"],
        "opdet_patterns":  [r"^volArapuca"],
        "struct_patterns": [r"^volTPCActive$", r"^volGaseousArgon$",
                            r"^volTPC$", r"^volCryostat$"],
        "hint_vols": {"volAnodePlate"},
    },
    "hd": {
        "label": "HD (Horizontal Drift / APA)",
        "anode_patterns":  [r"^volAPAFrameYSide$"],
        "opdet_patterns":  [r"^volArapuca_\d"],
        "struct_patterns": [r"^volTPCActiveInner$", r"^volGaseousArgon$",
                            r"^volTPC$", r"^volCryostat$"],
        "hint_vols": {"volAPAFrameYSide"},
    },
}

_GENERIC_PROFILE: dict = {
    "label": "DUNE (generic)",
    "anode_patterns":  [r"volAnode|volCRM$|APAFrame"],
    "opdet_patterns":  [r"volOpDetSensitive|volArapuca"],
    "struct_patterns": [r"volTPCActive$", r"volGaseousArgon$",
                        r"^volTPC$", r"^volCryostat$"],
}


def detect_profile(vol_names: set[str]) -> str | None:
    for key, prof in PROFILES.items():
        if prof["hint_vols"] & vol_names:
            return key
    return None


def get_profile(key: str | None) -> tuple[str, dict]:
    """Return (key, profile_dict); falls back to generic if key is None."""
    if key and key in PROFILES:
        return key, PROFILES[key]
    return "", _GENERIC_PROFILE
