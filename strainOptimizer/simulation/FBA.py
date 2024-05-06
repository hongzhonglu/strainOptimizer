# -*- coding: utf-8 -*-
# date : 2023/3/18 
# author : wangh
# file : FBA.py
# project : strainOptimizer

def fba(model, targetID, c_source, c_uptake,model_type='etfl',direction='max'):
    if model_type not in ['etfl', 'ecGEM','GEM']:
        raise ValueError('model_type should be either "etfl" or "ecGEM"')
    if c_uptake is None:
        c_uptake = 1
    if model_type == 'etfl':
        model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake
    elif model_type == 'ecGEM':
        model.reactions.get_by_id(c_source).bounds = c_uptake, c_uptake
    elif model_type == 'GEM':
        model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake

    model.objective = targetID
    model.objective_direction = direction
    sol = model.optimize()

    return sol