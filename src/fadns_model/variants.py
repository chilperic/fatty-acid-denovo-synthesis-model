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
