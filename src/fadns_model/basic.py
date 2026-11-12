"""Basic mass-action FAS model as a stoichiometric system.

Equivalent to the preserved ``basic_model.system``; no fitted parameter values
or physical units are inferred here. State order matches that legacy model.
"""

import numpy as np

COMPOUNDS = (
    "AcetylCoA",
    "MalonylCoA",
    *(f"C{n}FA" for n in range(4, 20, 2)),
    "E",
    *(f"EC{n}" for n in range(2, 20, 2)),
    "NADPH",
)
PARAMETERS = ("beta", "kappa", "delta", *(f"kd{i}" for i in range(1, 9)))
# Binding; eight elongations; eight chain-release reactions.
STOICHIOMETRY = np.zeros((21, 17))
STOICHIOMETRY[[0, 10, 11], 0] = [-1, -1, 1]
for j in range(8):
    STOICHIOMETRY[[1, 11 + j, 12 + j, 20], 1 + j] = [-1, -1, 1, -2]
    STOICHIOMETRY[[2 + j, 10, 12 + j], 9 + j] = [1, 1, -1]
STOICHIOMETRY.flags.writeable = False


def validate_parameters(params):
    if set(params) != set(PARAMETERS):
        raise ValueError(f"Expected exactly these basic-model parameters: {PARAMETERS}")
    if any(not np.isscalar(v) or not np.isfinite(v) or v < 0 for v in params.values()):
        raise ValueError("Mass-action constants must be finite and nonnegative")


def rhs(t, y, params):
    y = np.asarray(y, dtype=float)
    if y.shape != (21,):
        raise ValueError("Basic model requires 21 states in COMPOUNDS order")
    fluxes = np.empty(17)
    fluxes[0] = params["beta"] * y[10] * y[0]
    fluxes[1:9] = (
        y[11:19] * y[1] * y[20] * np.array([params["kappa"]] * 7 + [params["delta"]])
    )
    fluxes[9:] = y[12:20] * np.array([params[f"kd{i}"] for i in range(1, 9)])
    return STOICHIOMETRY @ fluxes
