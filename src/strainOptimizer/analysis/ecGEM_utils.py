# -*- coding: utf-8 -*-
# date : 2024/3/20 
# author : wangh
import pandas as pd


def prepare_prot_solution_for_ec(solution, enzymeIDlist=None):
    '''Extract the protein abundances result from solution for ecGEM model.
    parameters:
        solution: Cobra solution
        enzymeIDlist: a list of enzyme ID(optional)
    return:
        prots_solution: pd.Series, the protein abundances result
        '''
    prots_solution = pd.Series()

    if enzymeIDlist is None:
        for id in solution.fluxes.index:
            if 'draw_prot_' in id:
                prots_solution[id] = solution.fluxes[id]

    else:
        for enz in enzymeIDlist:
            prots_solution[enz] = solution.fluxes[enz]

    return prots_solution


def prepare_metabolic_solution_for_ec(solution, rxnList=None):
    '''Ectract all metabolic fluxes data result from solution for etfl model.
    parameters:
        solution: Cobra solution
        rxnList: a list of reaction ID
    return:
        fluxes: pd.Series, the fluxes data
    '''
    metabolic_solution = pd.Series()
    if rxnList is None:
        for id in solution.fluxes.index:
            if id.startswith('r_'):
                metabolic_solution[id] = solution.fluxes[id]
    else:
        for rxn in rxnList:
            if rxn in solution.fluxes.index:
                metabolic_solution[rxn] = solution.fluxes[rxn]

    return metabolic_solution