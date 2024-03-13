# -*- coding: utf-8 -*-
# date : 2024/3/12 
# author : wangh
from optlang.symbolics import Zero, add
from cobra.util import solver as sutil
from strainOptimizer.simulation.pprotFBA import pprotFBA_prot_conc
import pandas as pd
from pytfa.optim.utils import symbol_sum
import os


def mopa(model, reference_prots, linear: bool = True):
    """Compute a single solution based on (linear) Minimization of proteomic adjustment(MAPA) for a provided set of
    proteins.
    parameters:
        model: cobra.Model
        reference_prots: pd.Series, the reference protein abundances. This can only provide a customized set of proteins
         to be constrained.
        linear: bool, whether to use linear or quadratic optimization. Default is True.
    return:
        solution: the solution of MOPA
    """
    with model:
        add_mopa(model=model, reference_prots=reference_prots, linear=linear)
        print(model.objective.expression)
        solution = model.optimize()
    return solution


def add_mopa(model, reference_prots, linear: bool = True):
    """
    Add MOPA constraints and objective representing to the model.
    """
    if "mopa_old_objective" in model.solver.variables:
        raise ValueError("The model is already adjusted for MOPA.")

    # Fall back to default QP solver if current one has no QP capability
    if not linear and sutil.interface_to_str(model.problem) not in sutil.qp_solvers:
        model.solver = sutil.choose_solver(model, qp=True)

    prob = model.problem
    v = prob.Variable("mopa_old_objective")
    c = prob.Constraint(
        model.solver.objective.expression - v,
        lb=0.0,
        ub=0.0,
        name="mopa_old_objective_constraint",
    )
    to_add = [v, c]
    model.objective = prob.Objective(Zero, direction="min", sloppy=True)
    obj_vars = []
    all_enzymeList=reference_prots[reference_prots.index.str.contains('EZ_')].index
    for protID in all_enzymeList:
        enzID=protID.replace('EZ_','')
        prot = model.enzymes.get_by_id(enzID)
        ref_abundance = reference_prots[protID]
        if linear:
            components = sutil.add_absolute_expression(
                model,
                prot.variable,
                name="mopa_dist_" + prot.id,
                difference=ref_abundance,
                add=False,
            )
            to_add.extend(components)
            obj_vars.append(components.variable)
        else:
            dist = prob.Variable("mopa_dist_" + prot.id)
            # option 1: use scaled protein abundance(0<abundance<1)
            # const = prob.Constraint(
            #     prot.variable - dist,
            #     lb=ref_abundance,
            #     ub=ref_abundance,
            #     name="mopa_constraint_" + prot.id,
            # )
            # option 2: use the absolute protein abundance
            ref_abundance=ref_abundance*prot.scaling_factor
            const = prob.Constraint(
                prot.concentration - dist,
                lb=ref_abundance,
                ub=ref_abundance,
                name="mopa_constraint_" + prot.id,
            )
            to_add.extend([dist, const])
            obj_vars.append(dist ** 2)
    model.add_cons_vars(to_add)
    if linear:
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
    else:
        model.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)


def compare_mutant_vs_reference(ref, mutant, productID,threshold=5):
    """
    Compare the protein abundace change of mutant vs reference
    """
    to_ignore=['obj','mu','available_substrate','uptake']
    prot_change = pd.Series()
    flux_change= pd.Series()
    for r in ref.index:
        if r in to_ignore:
            continue
        elif 'EZ_' in r:
            prot_ref= ref[r]
            prot_mutant= mutant.raw[r]
            if prot_ref > 0:
                fc = prot_mutant / prot_ref
                prot_change[r] = fc
            elif prot_ref == 0:
                if prot_mutant > 0:
                    fc = 1000
                elif prot_mutant == 0:
                    fc = 1
                prot_change[r] = fc
        else:
        # if 'draw_prot_' in r:
            flux_ref = ref[r]
            flux_mutant = mutant[r]
            if flux_ref != 0:
                fc = flux_mutant / flux_ref
                flux_change[r] = fc
            elif flux_ref == 0:
                if flux_mutant != 0:
                    fc = 1000
                elif flux_mutant == 0:
                    fc = 1
                flux_change[r] = fc
    ref_product=ref[productID]
    mutant_product=mutant.raw[productID]
    print(f'product {productID} wt vs mutant: {ref_product} vs {mutant_product}')
    print('Protein abundace change more than 5 times:', len(prot_change[prot_change > threshold]))
    print('Protein abundace change less than 0.2 times:', len(prot_change[prot_change < 1 / threshold]))
    print('Flux change more than 5 times:', len(flux_change[flux_change > threshold]))
    print('Flux change less than 0.2 times:', len(flux_change[flux_change < 1 / threshold]))
    return prot_change, flux_change


def constrain_enzymes_based_abs_abundance(model, select_enzyme='EZ_TPI', ub0=1):
    from etfl.optim.variables import EnzymeVariable
    from etfl.optim.constraints import ModelConstraint
    # a function to add a constraint on total amount of enzymes based on their
    # fraction from total amount of proteins (should be before adding dummy)
    enz_vars = model.get_variables_of_type(EnzymeVariable)
    #enz_vars = [x for x in enz_vars if x.name not in exclusion]
    single_vars = [x for x in enz_vars if x.name == select_enzyme]
    expr = symbol_sum([x for x in single_vars])
    model.add_constraint(kind=ModelConstraint,
                         hook=model,
                         expr=expr,
                         id_='enzyme_fix_' + select_enzyme,
                         lb=ub0,  # cannot be negative
                         ub=ub0*1.1
                         ) # once this value is changed, it can't be replaced by the new value
    model.repair()
    return model


def prep_sol(substrate_uptake, model,GLC_RXN_ID ='EX_glc__D_e'):
    ret = {'obj':model.solution.objective_value,
           'mu':model.solution.fluxes.loc[model.growth_reaction.id],
           'available_substrate':-1*substrate_uptake,
           'uptake':-1*model.solution.fluxes[GLC_RXN_ID]
           }

#    for exch in model.exchanges:
#        ret[exch.id] = model.solution.fluxes.loc[exch.id]
    for rxn in model.reactions:
        ret[rxn.id] = model.solution.fluxes.loc[rxn.id]
    for enz in model.enzymes:
        ret['EZ_'+ enz.id] = model.solution.raw.loc['EZ_'+enz.id]
        ret=pd.Series(ret)
    return ret


if __name__=="__main__":
    from etfl.io.json import load_json_model, save_json_model

    solver = 'optlang-gurobi'

    ecoli = load_json_model("examples/models/ecoli/ecoli_core_curated.json", solver=solver)

    # set tolerance as 1e-9
    ecoli.tolerance = 1e-9

    # turn off by-product
    # ecoli.reactions.get_by_id('EX_for_e').bounds = 0, 0
    # ecoli.reactions.get_by_id('EX_ac_e').bounds = 0, 0
    # ecoli.reactions.get_by_id('EX_pyr_e').bounds = 0, 0
    # ecoli.reactions.get_by_id('EX_lac__D_e').bounds = 0, 0
    # ecoli.reactions.get_by_id('EX_etoh_e').bounds = 0, 0
    # ecoli.reactions.get_by_id('EX_acald_e').bounds = 0, 0

    model=ecoli
    productID='EX_succ_e'
    growth_id='Biomass_Ecoli_core'
    glucose_id='EX_glc__D_e'
    dummy_enzID='dummy_enzyme'

    # 1. get reference solution and enzyme concentrations
    model.objective = growth_id
    model.objective_direction = 'max'
    growth = model.slim_optimize()
    # fix growth and minimize total enzyme usage
    model.reactions.get_by_id(growth_id).bounds = growth, growth
    model.objective = model.enzymes.get_by_id(dummy_enzID).variable
    ref_solution = model.optimize()
    reference=prep_sol(substrate_uptake=-10, model=model)


    mutant=model.copy()
    # release all constraints
    mutant.reactions.get_by_id(growth_id).bounds = 0,1000
    # mutant 1: directly increase the production
    # with model:
    #     model.reactions.get_by_id(growth_id).bounds = 0,1000
    #     model.objective = productID
    #     max_product = model.slim_optimize()
    #     product=0.5*max_product
    # mutant.reactions.get_by_id(productID).bounds = product,product

    # mutant 2: overexpress a protein
    protID='PFK'
    # protID='G6PDH2r'
    # protID='PGI'
    test_pro_abundance = reference['EZ_' + protID]
    test_pro_abundance=test_pro_abundance.round(6)
    mutant = constrain_enzymes_based_abs_abundance(model=model, select_enzyme='EZ_' + protID,
                                                  ub0=test_pro_abundance * 5)
    mutant.objective = growth_id
    mutant.objective_direction = 'max'
    mutant.slim_optimize()

    # enz_constraint=mutant.get_constraints_of_type('ModelConstraint').get_by_id('enzyme_fix_EZ_'+protID)
    # mutant.remove_constraint(enz_constraint)


    mutant.objective=model.enzymes.get_by_id(protID).variable
    with mutant:
        mutant.objective_direction='max'
        print(mutant.slim_optimize())
    with mutant:
        mutant.objective_direction='min'
        print(mutant.slim_optimize())

    # test MOMA
    sol_linear = mopa(model=ecoli, reference_prots=reference, linear=True)
    sol= mopa(model=mutant, reference_prots=reference, linear=False)

    _,_=compare_mutant_vs_reference(reference, sol_linear, productID,threshold=5)
    _,_=compare_mutant_vs_reference(reference, sol, productID,threshold=5)
