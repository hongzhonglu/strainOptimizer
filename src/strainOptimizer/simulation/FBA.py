# -*- coding: utf-8 -*-
from .utils import set_carbon_source_bounds


def fba(model, target_id, c_source, c_uptake=None, model_type='etfl', direction='max'):
    if model_type not in ['etfl', 'ecGEM', 'GEM']:
        raise ValueError('model_type should be either "etfl", "ecGEM", or "GEM"')
    if c_uptake is None:
        c_uptake = 1
    set_carbon_source_bounds(model, c_source, c_uptake, model_type, fixed=True)
    model.objective = target_id
    model.objective_direction = direction
    sol = model.optimize()
    return sol
