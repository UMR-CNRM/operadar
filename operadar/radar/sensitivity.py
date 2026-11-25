#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Radar sensitivity: the simulated variables are removed where the reflectivity is below the minimum
detectable reflectivity of the radar.

The minimum detectable reflectivity increases with the range r (the received power decreases as 1/r^2):
    Zmin(r) = Zmin_1km + 20 log10(r / 1 km)    [dBZ]
Where Zh < Zmin(r), Zh, Zdr, Kdp and Rhohv (and Zh_att) are set to NaN: a real radar would not measure
them, and the polarimetric variables computed from an extremely small signal (rhohv especially) are
meaningless.

The range is known when the radar geometry is simulated (radarloc), and for a vertically pointing radar
(cnst_angle = 90: |altitude - radar_altitude|). Otherwise (horizontal beam at every grid point, no
range) the threshold is Zmin_1km everywhere.
"""

import numpy as np

# Default minimum detectable reflectivity at 1 km (dBZ) per band. These are orders of magnitude for
# typical radars (operational L to X band radars, Ku to Ka precipitation radars, W band cloud radars):
# set Zmin_1km in the configuration file for a given radar.
DEFAULT_ZMIN_1KM = {"L": -30.0, "S": -30.0, "C": -30.0, "X": -30.0,
                    "Ku": -40.0, "K": -40.0, "Ka": -40.0, "W": -50.0}

MASKED_VARIABLES = ["Zh", "Zh_att", "Zdr", "Kdp", "Rhohv"]


def minimum_detectable_reflectivity(zmin_1km: float, slant_range: np.ndarray | None,
                                    shape: tuple) -> np.ndarray:
    """Zmin (dBZ) at each grid point; constant Zmin_1km where the range is unknown (None or NaN)."""
    zmin = np.full(shape, zmin_1km, dtype=float)
    if slant_range is not None:
        r = np.broadcast_to(slant_range, shape)
        known = np.isfinite(r)
        zmin[known] = zmin_1km + 20.0 * np.log10(np.maximum(r[known], 1.0) / 1000.0)
    return zmin


def vertical_range(altitude: np.ndarray, radar_altitude: float, shape: tuple) -> np.ndarray:
    """Range of a vertically pointing radar (zenith or nadir) at each grid point (m)."""
    return np.abs(np.broadcast_to(altitude, shape) - radar_altitude)


def apply_sensitivity_threshold(dpolDict: dict[str, np.ndarray], band: str, zmin_1km: float | None = None,
                                slant_range: np.ndarray | None = None, verbose: bool = True,
                                ) -> dict[str, np.ndarray]:
    """Set Zh, Zh_att, Zdr, Kdp and Rhohv to NaN where the reflectivity seen by the radar (Zh_att if the
    attenuation was computed, else Zh) is below the minimum detectable reflectivity."""
    zkey = "Zh_att" if "Zh_att" in dpolDict else "Zh"
    if zkey not in dpolDict:
        print("\tSensitivity threshold not applied: Zh is not computed (add 'Zh' to dpol2add)")
        return dpolDict
    if zmin_1km is None:
        zmin_1km = DEFAULT_ZMIN_1KM[band]
    zh = dpolDict[zkey]
    zmin = minimum_detectable_reflectivity(zmin_1km, slant_range, zh.shape)
    below = ~(zh >= zmin)  # also True where Zh is NaN
    for var in MASKED_VARIABLES:
        if var in dpolDict:
            dpolDict[var] = np.where(below, np.nan, dpolDict[var])
    if verbose:
        detected = np.isfinite(zh) & ~below
        print(f"\tSensitivity threshold: Zmin = {zmin_1km:g} dBZ at 1 km"
              + (" + 20 log10(r / 1 km)" if slant_range is not None else " (no range: constant)")
              + f", {int(np.sum(np.isfinite(zh) & below))} points removed, {int(np.sum(detected))} kept")
    return dpolDict
