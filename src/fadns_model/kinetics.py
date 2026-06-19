"""Generic kinetic helper functions shared by model variants."""

from __future__ import annotations

import numpy as np


def mass_action(k: float, *substrates: float) -> float:
    rate = float(k)
    for s in substrates:
        rate *= float(s)
    return rate


def reversible_mass_action(kf: float, kb: float, substrates, products) -> float:
    forward = float(kf)
    backward = float(kb)
    for s in substrates:
        forward *= float(s)
    for p in products:
        backward *= float(p)
    return forward - backward


def michaelis_menten(vmax: float, km: float, substrate: float) -> float:
    substrate = float(substrate)
    return float(vmax) * substrate / (float(km) + substrate)


def conservation(total: float, *values: float) -> float:
    return float(total) - float(np.sum(values))
