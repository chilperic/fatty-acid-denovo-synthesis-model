VARIANTS = {
    "coa_inhibition": "CoA-inhibition",
    "coa_inhibition_only_14_16_18_can_be_produced": "CoA-inhibition_only_14_16_18_can_be_produced",
    "evaluation": "Evaluation",
    "model_with_coa_inhibition_rate_of_release_from_distribution": "Model_with_CoA-inhibition_rate_of release_from _distribution",
    "basic_model": "basic-model",
    "basic_model_only_14_16_18": "basic-model_only_14_16_18",
    "malonyl_coa_and_acetyl_coa_competition": "malonyl-coa-and-acetyl-coa-competition",
    "malonyl_coa_competition": "malonyl-coa-competition",
    "open_system_coa_inhibition_only_14_16_18_can_be_produced": "open_system_CoA-inhibition_only_14_16_18_can_be_produced",
    "sinking_enzyme": "sinking-enzyme",
}


def list_variants() -> dict[str, str]:
    return dict(VARIANTS)


def load_variant(name):
    """Load a preserved equation implementation by its explicit registry name.

    Legacy plotting/fitting helpers are not part of the supported modern API.
    Use ``simulate_variant`` instead of the legacy Assimulo integration wrapper.
    """
    from importlib import import_module

    if name not in VARIANTS:
        raise ValueError(f"Unknown variant {name!r}; choose from {list(VARIANTS)}")
    return import_module(f".legacy_variants.{name}", package="fadns_model")


def simulate_variant(name, y0, params, **solver_options):
    """Run the original variant equations with validated SciPy integration."""
    import numpy as np

    from .simulation import simulate_time_course

    module = load_variant(name)
    if np.asarray(y0).shape != (len(module.compounds),):
        raise ValueError("y0 must follow the variant's compounds order")
    result = simulate_time_course(module.system, y0, params, **solver_options)
    result.compounds = tuple(module.compounds)
    result.variant = name
    return result
