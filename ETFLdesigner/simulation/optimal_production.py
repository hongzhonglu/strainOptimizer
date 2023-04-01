# -*- coding: utf-8 -*-
# date : 2023/3/23 
# author : wangh

from etfl.optim.utils import safe_optim
from ETFLdesigner.ETFLdesigner.simulation import pprotFBA


def optim_production_simulating(model, target, c_source,c_uptake=1, alpha=1, tol=1e-6,model_type='etfl'):
    """
    Function to simulate flux distribution in the optimal production condition.
    Performing a series of LP optimizations on an ETFL/ecModel,
    by first maximizing biomass, then fixing a suboptimal value and
    proceeding to maximize a given production target reaction, last
    a protein pool minimization is performed subject to the optimal
    production level.

    Args:
    - model: (cobra.Model) the model to optimize
    - target: (str) reaction ID for the objective reaction
    - c_source: (str) reaction ID for the main carbon source uptake reaction
    - c_uptake: (float) uptake rate for the main carbon source (default: 1)
    - alpha: (float) scaling factor for desired suboptimal growth (default: 1)
    - tol: (float) numerical tolerance for fixing bounds (default: 1e-6)
    - model_type: (str) type of model to optimize.'etfl'/'ecGEM' (default: 'etfl')

    Returns:
    - flux: (pandas.Series) a Series of the optimized flux distribution

    """
    # Fix a unit main carbon source uptake
    # model=model.copy()
    uptake_rxn = model.reactions.get_by_id(c_source)
    # c_source_MW = uptake_rxn.reactants[0].formula_weight/1000
    if model_type=='etfl':
        uptake_rxn.bounds = -c_uptake, -c_uptake
        # Maximize growth
        gr_rxn=model.growth_reaction
        grID=gr_rxn.id
        model.objective = grID
        sol = safe_optim(model)
        max_gr=sol.fluxes[grID]

    elif model_type=='ecGEM':
        uptake_rxn.bounds = c_uptake - tol, c_uptake + tol
        # Maximize growth
        gr_rxn = model.reactions.get_by_id('r_2111')
        grID = gr_rxn.id
        model.objective = grID
        sol=model.optimize()
        max_gr = sol.fluxes[grID]

    # Fix growth suboptimal and then maximize product
    model.reactions.get_by_id(grID).lower_bound = max_gr * (1 - tol) * alpha
    flux = pprotFBA.ppFBA(model, target,c_source=c_source,c_uptake=c_uptake, tol=tol,model_type=model_type)


    # restore the original bounds
    model.reactions.get_by_id(grID).lower_bound = 0

    return flux
