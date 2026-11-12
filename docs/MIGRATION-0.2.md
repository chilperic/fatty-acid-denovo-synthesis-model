# Scientific and API migration to 0.2.0

The modern basic model is expressed as a stoichiometric matrix and mass-action flux vector. Random-state derivative tests verify equivalence to the original basic equations; enzyme and acetyl-unit conservation are tested. Other named variants use an explicit registry and validated SciPy integration, retaining historical equations. Legacy plotting and fitting helpers contain exploratory code and are outside the supported API.

CSV loading now retains time zero in headerless files. Small negative digitized observations are preserved. Headerless time units and all measurement units require source verification; do not silently infer concentration units. Synthetic example parameters are not a paper fit. Reproduction of every published figure and all legacy fitted parameter sets remains pending.
