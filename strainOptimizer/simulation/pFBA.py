# -*- coding: utf-8 -*-
from cobra.flux_analysis import pfba
from strainOptimizer.analysis import prepare_metabolic_solution_for_etfl, prepare_metabolic_solution_for_ec

def pFBA(model, targetID, c_source, c_uptake, model_type='etfl', direction='max'):
    if model_type not in ['etfl', 'ecGEM', 'GEM']:
        raise ValueError('model_type should be either "etfl" or "ecGEM"')
    if c_uptake is None:
        c_uptake = 1
    if model_type == 'etfl':
        model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake
        met_rxnList=[rxn for rxn in model.reactions if rxn.id.startswith('r_')]
    elif model_type == 'ecGEM':
        model.reactions.get_by_id(c_source).bounds = c_uptake, c_uptake
        met_rxnList=[rxn for rxn in model.reactions if rxn.id.startswith('r_')]
    elif model_type == 'GEM':
        model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake
        met_rxnList=[rxn for rxn in model.reactions]


    model.objective = targetID
    model.objective_direction = direction
    sol = pfba(model=model,
               fraction_of_optimum=1.0,
               reactions=met_rxnList,
               )

    return sol