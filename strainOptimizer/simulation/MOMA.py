# -*- coding: utf-8 -*-
# date : 2024/3/20 
# author : wangh
'''Minimization of the metabolic adjustment(MOMA) accoding to the reference rxneins resourec allocation.
'''
from optlang.symbolics import Zero, add
from cobra.util import solver as sutil
import pandas as pd
from strainOptimizer.analysis import prepare_metabolic_solution_for_etfl, prepare_metabolic_solution_for_ec


def moma(model, reference_solution, linear: bool = True,model_type='ecGEM'):
    """Compute a solution based on (linear) Minimization of metabolic adjustment(MOMA) for a provided set of rxns.
    parameters:
        model: ETFL model or ecGEM model
        reference_solution: reference_solution: reference optimization solution by cobrapy or pyTFA.
        linear: bool, whether to use linear or quadratic optimization. Default is True.
        model_type: str, 'ecGEM' or 'etfl'
    return:
        solution: the solution of MOMA
    """
    with model:
        if model_type=='etfl':
            reference_fluxes= prepare_metabolic_solution_for_etfl(reference_solution)
            add_moma_etfl(model=model, reference_fluxes=reference_fluxes, linear=linear)
        elif model_type=='ecGEM':
            reference_fluxes= prepare_metabolic_solution_for_ec(reference_solution)
            add_moma_ecGEM(model=model, reference_fluxes=reference_fluxes, linear=linear)
        else:
            raise ValueError('model_type should be "etfl" or "ecGEM"')
        # print(model.objective.expression)
        solution = model.optimize()
    return solution



def add_moma_etfl(model, reference_fluxes, linear: bool = True):
    """
    Add MOMA constraints and objective representing to the ETFL model.
    parameters:
        model: ETFL model
        reference_fluxes: pd.Series, the reference rxnein abundances. This can only provide a customized set of rxneins
         to be constrained.
        linear: bool, whether to use linear or quadratic optimization. Default is True.
    """

    if "moma_old_objective" in model.solver.variables:
        raise ValueError("The model is already adjusted for MOMA.")

    # Fall back to default QP solver if current one has no QP capability
    if not linear and sutil.interface_to_str(model.problem) not in sutil.qp_solvers:
        model.solver = sutil.choose_solver(model, qp=True)

    prob = model.problem
    v = prob.Variable("moma_old_objective")
    c = prob.Constraint(
        model.solver.objective.expression - v,
        lb=0.0,
        ub=0.0,
        name="moma_old_objective_constraint",
    )
    to_add = [v, c]
    model.objective = prob.Objective(Zero, direction="min", sloppy=True)
    obj_vars = []
    for rxnID in reference_fluxes.index:
        rxn=model.reactions.get_by_id(rxnID)
        ref_flux = reference_fluxes[rxnID]
        if linear:
            components = sutil.add_absolute_expression(
                model,
                rxn.flux_expression,
                name="moma_dist_" + rxnID,
                difference=ref_flux,
                add=False,
            )
            to_add.extend(components)
            obj_vars.append(components.variable)
        else:
            dist = prob.Variable("moma_dist_" + rxn.id)
            const = prob.Constraint(
                rxn.flux_expression - dist,
                lb=ref_flux,
                ub=ref_flux,
                name="moma_constraint_" + rxn.id,
            )
            to_add.extend([dist, const])
            obj_vars.append(dist ** 2)
    model.add_cons_vars(to_add)

    # set MOMA objective
    if linear:
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
    else:
        model.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)


def add_moma_ecGEM(model, reference_fluxes, linear: bool = True):
    """
    Add MOMA constraints and objective representing to the ecGEM model.
    """
    if "moma_old_objective" in model.solver.variables:
        raise ValueError("The model is already adjusted for MOMA.")

    # Fall back to default QP solver if current one has no QP capability
    if not linear and sutil.interface_to_str(model.problem) not in sutil.qp_solvers:
        model.solver = sutil.choose_solver(model, qp=True)

    prob = model.problem
    v = prob.Variable("moma_old_objective")
    c = prob.Constraint(
        model.solver.objective.expression - v,
        lb=0.0,
        ub=0.0,
        name="moma_old_objective_constraint",
    )
    to_add = [v, c]
    model.objective = prob.Objective(Zero, direction="min", sloppy=True)
    obj_vars = []
    for rxnID in reference_fluxes.index:
        rxn=model.reactions.get_by_id(rxnID)
        ref_abundance = reference_fluxes[rxnID]
        if linear:
            components = sutil.add_absolute_expression(
                model,
                rxn.flux_expression,
                name="moma_dist_" + rxn.id,
                difference=ref_abundance,
                add=False,
            )
            to_add.extend(components)
            obj_vars.append(components.variable)
        else:
            dist = prob.Variable("moma_dist_" + rxn.id)
            const = prob.Constraint(
                rxn.flux_expression - dist,
                lb=ref_abundance,
                ub=ref_abundance,
                name="moma_constraint_" + rxn.id,
            )
            to_add.extend([dist, const])
            obj_vars.append(dist**2)
    model.add_cons_vars(to_add)
    if linear:
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
    else:
        model.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)
