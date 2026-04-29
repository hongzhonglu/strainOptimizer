# -*- coding: utf-8 -*-
# date : 2024/3/20 
# author : wangh
import pandas as pd


def prepare_prot_solution_for_etfl(solution,enzymeIDlist=None):
    '''Extract the protein abundances result from solution for ecGEM model.
    parameters:
        solution: Cobra solution
        enzymeIDlist: a list of enzyme ID(optional)
    return:
        prots_solution: pd.Series, the protein abundances result
        '''
    prots_solution = pd.Series()

    if enzymeIDlist is None:
        for id in solution.raw.index:
            if id.startswith('EZ_'):
                prots_solution[id]=solution.raw[id]

    else:
        for enz in enzymeIDlist:
            prots_solution[enz]=solution.raw[enz]

    return prots_solution


def prepare_metabolic_solution_for_etfl(solution, rxnList=None, flux_tol=1e-6):
    '''Extract metabolic fluxes from an ETFL solution for use as MOMA reference.

    Only reactions with |flux| > flux_tol are returned. Filtering out near-zero
    fluxes avoids over-constraining the MOMA problem when the perturbed model
    has a slightly different feasible region.

    parameters:
        solution: pyTFA solution
        rxnList: a list of reaction IDs (None = all r_ reactions)
        flux_tol: minimum absolute flux to include (default 1e-6)
    return:
        fluxes: pd.Series, the fluxes data
    '''
    metabolic_solution = pd.Series(dtype=float)
    if rxnList is None:
        for id in solution.fluxes.index:
            if id.startswith('r_') and abs(solution.fluxes[id]) > flux_tol:
                metabolic_solution[id] = solution.fluxes[id]
    else:
        for rxn in rxnList:
            if rxn in solution.fluxes.index and abs(solution.fluxes[rxn]) > flux_tol:
                metabolic_solution[rxn] = solution.fluxes[rxn]

    return metabolic_solution

