# -*- coding: utf-8 -*-
# date : 2023/3/18 
# author : wangh
import numpy as np
import pandas as pd
from etfl.io.json import load_json_model
from tqdm import tqdm
from pytfa.analysis.variability import _variability_analysis_element
from etfl.optim.utils import safe_optim
from cobra.flux_analysis import flux_variability_analysis


def etfl_EVA(model,targetID,enzymeIDlist,c_source,c_uptake,fraction_of_optimum=0.99,obj_direction='max'):
    '''do enzyme variety analysis for ETFL model
    para:
        model: ETFL model
        enzymeIDlist: a list of enzyme ID
        fraction_of_optimum: Requires that the objective value is at least the
            fraction times maximum objective value.Must be <= 1.0. (default 0.95)
    return:
        a dataframe of FVA result
        '''
    # fix substrate uptake
    model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake

    # 1.Optimize a given objective
    model.objective = targetID
    model.objective_direction = obj_direction
    sol = safe_optim(model)
    obj_value = sol.objective_value

    # 2. get all enzyme variable
    all_enz = model.get_variables_of_type('EnzymeVariable')
    all_enzIDlist = [enz.id for enz in all_enz]
    # get the target enzyme list
    target_enzlist = {}
    for enzID in enzymeIDlist:
        if enzID in all_enzIDlist:
            target_enzlist[enzID] = model.enzymes.get_by_id(enzID).variable
        else:
            print(f"can't find Enzyme {enzID} in the {model.name}")

    # 3. fix old objective value and add constraint
    if model.solver.objective.direction == "max":
        fva_old_objective = model.problem.Variable(
            "fva_old_objective",
            lb=fraction_of_optimum * obj_value,
        )
    else:
        fva_old_objective = model.problem.Variable(
            "fva_old_objective",
            ub=fraction_of_optimum * obj_value,
        )
    fva_old_obj_constraint = model.problem.Constraint(
        model.solver.objective.expression - fva_old_objective,
        lb=0,
        ub=0,
        name="fva_old_objective_constraint",
    )
    model.add_cons_vars([fva_old_objective, fva_old_obj_constraint])

    # 5.do enzyme variety analysis
    results = {'min': {}, 'max': {}}
    for sense in ['min', 'max']:
        for k, var in tqdm(target_enzlist.items(), desc=sense + 'imizing'):
            model.logger.debug(sense + '-' + k)
            results[sense][k] = _variability_analysis_element(model, var, sense)

    # 6.remove fixed constraint and old objective
    model.remove_cons_vars([fva_old_objective, fva_old_obj_constraint])
    # restore old objective
    model.objective = targetID

    df = pd.DataFrame(results)
    df.rename(columns={'min': 'minimum', 'max': 'maximum'}, inplace=True)

    return df


def ecGEM_EVA(model,targetID,enzymeIDlist,c_source,c_uptake,fraction_of_optimum=1,obj_direction='max'):
    '''do enzyme variety analysis for ecGEM
       para:
           model: ecGEM model
           enzymeIDlist: a list of enzyme ID
           fraction_of_optimum: Requires that the objective value is at least the
               fraction times maximum objective value.Must be <= 1.0. (default 0.95)
       return:
           a dataframe of FVA result
           '''
    # fix substrate uptake
    model.reactions.get_by_id(c_source).bounds = c_uptake, c_uptake

    # set the objective function
    model.objective = targetID
    model.objective_direction = obj_direction

    df_fva_result=flux_variability_analysis(model=model,reaction_list=enzymeIDlist,fraction_of_optimum=fraction_of_optimum)

    return df_fva_result


def enzymeVA(model,targetID,enzymeIDlist,c_source,c_uptake,fraction_of_optimum=0.99,obj_direction='max',model_type='etfl'):
    '''do enzyme variety analysis for ecGEM/ETFL model
    para:
        model: ecGEM/ETFL model
        targetID: the target reaction ID
        enzymeIDlist: a list of enzyme ID
        c_source: the carbon source ID
        c_uptake: the carbon source uptake rate(default=1 mmol/gDW/h)
        fraction_of_optimum: Requires that the objective value is at least the
            fraction times maximum objective value.Must be <= 1.0. (default 0.99)
        obj_direction: the direction of the objective function(default='max')
        model_type: the type of the model(default='etfl')

    return:
        a dataframe of FVA result
        '''

    if model_type=='etfl':
        eva_result= etfl_EVA(model=model,targetID=targetID,enzymeIDlist=enzymeIDlist,c_source=c_source,c_uptake=c_uptake,fraction_of_optimum=fraction_of_optimum,obj_direction=obj_direction)
    elif model_type=='ecGEM':
        eva_result=ecGEM_EVA(model=model,targetID=targetID,enzymeIDlist=enzymeIDlist,c_source=c_source,c_uptake=c_uptake,fraction_of_optimum=fraction_of_optimum,obj_direction=obj_direction)

    return eva_result