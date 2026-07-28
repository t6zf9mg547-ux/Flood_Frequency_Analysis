"""
Bulletin 17B-style weighted skew for Pearson III / Log-Pearson III.

Station skew (computed from a single record) is notoriously unstable with
typical record lengths (n ~ 30-100 years) -- Bulletin 17B (Interagency
Advisory Committee on Water Data, 1982) addresses this by blending the
station skew with an independently-derived *regional* skew estimate,
weighting each by the inverse of its mean square error (MSE): a more
precise estimate gets more weight.

This module implements the classical Bulletin 17B procedure (station
skew + regional skew, weighted by MSE) -- NOT the full Bulletin 17C
Expected Moments Algorithm (EMA), which re-derives the weighted skew
iteratively alongside historical/censored-data handling. That's a
substantially larger undertaking; this covers the well-documented,
verifiable classical case.

Regional skew and its MSE are NOT something this tool can look up or
guess -- they come from a published regional study specific to your area
(e.g. a state DOT report, a USGS regional-skew study, or the Bulletin 17B
national skew map). You must supply both explicitly.

References:
  - Interagency Advisory Committee on Water Data (1982), Bulletin 17B,
    Appendix 8 (station skew MSE) and Equation 6 (weighted skew).
  - If using the Bulletin 17B national skew map for the regional value,
    its documented MSE is 0.302 (commonly superseded by more precise
    state/regional studies -- e.g. Texas DOT reports 0.123, Arizona USGS
    studies report as low as 0.08).
"""
from __future__ import annotations
import numpy as np


def station_skew_mse(skew: float, n: int) -> float:
    """
    Mean square error of the station (at-site) skew, as a function of its
    own magnitude and the record length -- Bulletin 17B Appendix 8,
    Equation 4.17(a-e) (as reproduced in numerous state DOT/USGS
    hydrology references):

        log10(MSE) = A - B * log10(n / 10)

        A = -0.33 + 0.08*|G|   if |G| <= 0.90
        A = -0.52 + 0.30*|G|   if |G| >  0.90
        B =  0.94 - 0.26*|G|   if |G| <= 1.50
        B =  0.55              if |G| >  1.50
    """
    G = abs(skew)
    A = (-0.33 + 0.08 * G) if G <= 0.90 else (-0.52 + 0.30 * G)
    B = (0.94 - 0.26 * G) if G <= 1.50 else 0.55
    log_mse = A - B * np.log10(n / 10.0)
    return float(10 ** log_mse)


def weighted_skew(station_skew: float, n: int, regional_skew: float,
                   regional_mse: float = 0.302) -> dict:
    """
    Combine station and regional skew per Bulletin 17B Equation 6:

        G_w = (MSE_regional * G_station + MSE_station * G_regional)
              / (MSE_station + MSE_regional)

    regional_mse defaults to 0.302, the documented MSE for the Bulletin
    17B national skew map -- but a region/state-specific study will
    almost always give a smaller (more precise) MSE and should be
    preferred if available.
    """
    mse_station = station_skew_mse(station_skew, n)
    Gw = ((regional_mse * station_skew + mse_station * regional_skew)
          / (mse_station + regional_mse))

    flag = None
    if abs(station_skew - regional_skew) > 0.5:
        flag = ("Station and regional skew differ by more than 0.5 -- Bulletin 17B "
                "recommends reviewing the data and basin characteristics in this case; "
                "consider whether the station skew (data quality/outliers) or the "
                "regional value (applicability to this basin) is more suspect.")

    return {
        "station_skew": float(station_skew), "station_mse": mse_station,
        "regional_skew": float(regional_skew), "regional_mse": float(regional_mse),
        "weighted_skew": float(Gw), "review_flag": flag,
    }
