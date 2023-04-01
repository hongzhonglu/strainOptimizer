# -*- coding: utf-8 -*-
# date : 2023/3/18 
# author : wangh
# file : pprotFBA.py
# project : ETFLdesigner

import os
import numpy as np
import pandas as pd
from pytfa.optim.utils import symbol_sum
from cobra.util.solver import set_objective
from etfl.optim.utils import safe_optim


def ppFBA(model, target,c_source,c_uptake=1,model_type='etfl',tol_ratio=0.01):
    """
   optimize for a given objective. Firstly,maximize the objective reaction,and then  a protein pool minimization is performed subject to the optimal
    production level.

    Args:
    - model: (cobra.Model) the model to optimize
    - target: (str) reaction ID for the objective reaction
    - c_source: (str) reaction ID for the main carbon source uptake reaction
    - c_uptake: (float) uptake rate for the main carbon source (default: 1)
    - tol: (float) numerical tolerance for fixing bounds(default: 1e-6)
    - prot_conc_output: (bool) whether to output protein concentration (default: False)
    - model_type: (str) type of model to optimize:'etfl'/'ecGEM' (default: 'etfl')

    Returns:
    - sol.x: (pandas.DataFrame) the optimized flux distribution

    """
    if model_type=='etfl':
        model.reactions.get_by_id(c_source).bounds = -c_uptake , -c_uptake
    elif model_type=='ecGEM':
        model.reactions.get_by_id(c_source).bounds = c_uptake, c_uptake
    # 1.set the target objective
    model.objective = target
    model.objective_direction = 'max'
    sol = safe_optim(model)
    max_obj = sol.objective_value

    # 2.Fix optimal value for objective and minimize total protein usage
    if model.solver.status == 'optimal':
        model.reactions.get_by_id(target).bounds = max_obj*(1-tol_ratio), max_obj
        # 3.minimize the total proteins to reuduce the solution space
        # for etfl model: maximize dummy enzyme pseudo rxn
        if model_type=='etfl':
            obj_expr = symbol_sum([model.enzymes.dummy_enzyme.variable])
            set_objective(model, obj_expr)
            model.objective_direction = 'max'
            sol2 = safe_optim(model)
        # for ecGEM model: minimize protein pool pseudo rxn
        elif model_type=='ecGEM':
            prot_pool_rxnID='prot_pool_exchange'
            model.objective = prot_pool_rxnID
            model.objective_direction = 'min'
            sol2=model.optimize()
    else:
        sol2 = np.zeros(len(model.reactions))
        return sol2

    # release the modified constraints and reset the objective as growth
    if model_type=='etfl':
        model.reactions.get_by_id(target).bounds = 0, 1000
        model.objective= model.growth_reaction.id
        model.objective_direction = 'max'
    elif model_type=='ecGEM':
        model.reactions.get_by_id(target).bounds = 0, 1000
        model.objective = 'r_2111'
        model.objective_direction = 'max'

    return sol2.fluxes


def etfl_ppFBA_prot_conc(model, target,c_source,c_uptake=1):
    '''use minprotFBA to predict target proteins concentration(notice!! the output is scaled protein concentration)
    para:
        model: must be ETFL model
        target: the target reaction ID
        enzID_list: a list of enzyme ID
        c_source: the carbon source ID
        c_uptake: the carbon source uptake rate(default=1 mmol/gDW/h)
        tol: the tolerance of the model
    return:
        a pandas series of protein concentration
        '''
    model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake
    # 1.Optimize for a given objective
    model.objective = target
    model.objective_direction = 'max'
    sol = safe_optim(model)
    max_obj = sol.objective_value
    all_enzIDlist= [enz.id for enz in model.enzymes]

    # 2.Fix optimal value for objective and minimize total protein usage
    if sol.status == 'optimal':
        model.reactions.get_by_id(target).bounds = max_obj, max_obj
        # for etfl model: minimize enzyme usage by maxing dummy enzyme
        obj_expr = symbol_sum([model.enzymes.dummy_enzyme.variable])
        set_objective(model, obj_expr)
        model.objective_direction = 'max'
        safe_optim(model)
        # get protein concentration
        prot_conc = pd.Series()
        for enzID in all_enzIDlist:
            prot_conc[enzID] = model.enzymes.get_by_id(enzID).scaled_X
        # restore the original bounds and objective
        model.reactions.get_by_id(target).bounds = 0, 1000
        model.objective = target
        model.objective_direction = 'max'
        return prot_conc

    else:
        # print the error message
        raise Exception('The model cannot be solved to optimality')


def ecGEM_ppFBA_prot_conc(model, target,c_source,c_uptake=1):
    '''use minprotFBA to predict target proteins concentration(1/3600 mmol/gDW)
    para:
        model: must be ecGEM model
        target: the target reaction ID
        enzID_list: a list of enzyme ID
        c_source: the carbon source ID
        c_uptake: the carbon source uptake rate(default=1 mmol/gDW/h)
        tol: the tolerance of the model
        return:
        a pandas series of protein concentration'''

    fluxes=ppFBA(model=model, target=target,c_source=c_source,c_uptake=c_uptake,model_type='ecGEM')
    all_prot_rxnIDlist=[rxn.id for rxn in model.reactions if rxn.id.startswith('draw_prot_')]
    all_prot_conc=fluxes[all_prot_rxnIDlist]

    return all_prot_conc

