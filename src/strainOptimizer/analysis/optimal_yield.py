# -*- coding: utf-8 -*-
from strainOptimizer.simulation.utils import set_carbon_source_bounds


def cal_max_yield(model, target_id, c_source, model_type='etfl', tol_ratio=0.001):
    '''calculate the maximum yield of the target product'''
    # 1. release carbon source uptake
    set_carbon_source_bounds(model, c_source, 1000, model_type, fixed=False)
    # 2. calculate maximum target product production
    model.objective = target_id
    model.objective_direction = 'max'
    max_prod = model.slim_optimize()
    if model.solver.status != 'optimal':
        return 0, 0
    if max_prod < 1e-9:
        return 0, 0
    # 3. fix the max target product production and minimize the C source uptake
    model.reactions.get_by_id(target_id).bounds = max_prod * (1 - tol_ratio), max_prod
    if model_type == 'etfl':
        model.objective = c_source
        model.objective_direction = 'max'
        opt_c_uptake = -model.slim_optimize()
    elif model_type == 'ecGEM':
        model.objective = c_source
        model.objective_direction = 'min'
        opt_c_uptake = model.slim_optimize()

    # 4. calculate the maximum yield
    try:
        max_yield = max_prod / opt_c_uptake
    except:
        max_yield = 0

    # restore the fixed bounds
    model.reactions.get_by_id(target_id).bounds = 0, 1000

    return max_yield, max_prod
