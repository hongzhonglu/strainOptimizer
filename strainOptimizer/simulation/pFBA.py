# -*- coding: utf-8 -*-
from cobra.flux_analysis import pfba
from strainOptimizer.analysis import prepare_metabolic_solution_for_etfl, prepare_metabolic_solution_for_ec

def pFBA(model, targetID, c_source, c_uptake, model_type='etfl', direction='max'):
    with model:
        if model_type not in ['etfl', 'ecGEM', 'GEM', 'GAN_ec']:
            raise ValueError('model_type should be either "etfl" or "ecGEM"')
        if c_uptake is None:
            c_uptake = 1
        if model_type == 'etfl':
            model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake
            met_rxnList=[rxn for rxn in model.reactions if rxn.id.startswith('r_')]
        elif model_type == 'ecGEM':
            model.reactions.get_by_id(c_source).bounds = c_uptake*(1-0.01), c_uptake*(1 + 0.01)  # allow a small tolerance for uptake
            # met_rxnList=[rxn for rxn in model.reactions if rxn.id.startswith('r_')]
            met_rxnList=[rxn for rxn in model.reactions if 'prot' not in rxn.id and not rxn.id.startswith('arm_')]
        elif model_type == 'GEM':
            model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake
            met_rxnList=[rxn for rxn in model.reactions]
        elif model_type == 'GAN_ec':
            model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake
            met_rxnList=[rxn for rxn in model.reactions if rxn.id.startswith('r_')]

        model.objective = targetID
        model.objective_direction = direction
        sol = pfba(model=model,
                   fraction_of_optimum=0.9,
                   reactions=met_rxnList,
                   )

    return sol