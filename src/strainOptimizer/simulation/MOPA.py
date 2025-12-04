# -*- coding: utf-8 -*-
# date : 2024/3/20 
# author : wangh
'''Minimization of the proteimic adjustment(MOPA) accoding to the reference proteins resourec allocation.
'''
from optlang.symbolics import Zero, add
from cobra.util import solver as sutil
import pandas as pd
from strainOptimizer.analysis import prepare_prot_solution_for_etfl, prepare_prot_solution_for_ec



def mopa(model, reference_solution, linear: bool = True,model_type='ecGEM',show=False):
    """Compute a solution based on (linear) Minimization of proteomic adjustment(MAPA) for a provided set of
    proteins.
    parameters:
        model: ETFL model or ecGEM model
        reference_solution: reference optimization solution by cobrapy or pyTFA.
        linear: bool, whether to use linear or quadratic optimization. Default is True.
        model_type: str, 'ecGEM' or 'etfl'
    return:
        solution: the solution of MOPA
    """
    model.tolerance = 1e-9
    with model:
        if model_type=='etfl':
            reference_prots= prepare_prot_solution_for_etfl(reference_solution)
            add_mopa_etfl(model=model, reference_prots=reference_prots, linear=linear)
        elif model_type=='ecGEM' or model_type=='GAN_ec':
            reference_prots= prepare_prot_solution_for_ec(reference_solution)
            add_mopa_ecGEM(model=model, reference_prots=reference_prots, linear=linear)
        else:
            raise ValueError('model_type should be "etfl" or "ecGEM"')
        if show:
            print(model.objective.expression)
        solution = model.optimize()
    return solution



def add_mopa_etfl(model, reference_prots, linear: bool = True):
    """
    Add MOPA constraints and objective representing to the ETFL model.
    parameters:
        model: ETFL model
        reference_prots: pd.Series, the reference protein abundances. This can only provide a customized set of proteins
         to be constrained.
        linear: bool, whether to use linear or quadratic optimization. Default is True.
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
    all_enzymeList = reference_prots[reference_prots.index.str.contains('EZ_')].index
    for protID in all_enzymeList:
        enzID = protID.replace('EZ_', '')
        prot = model.enzymes.get_by_id(enzID)
        ref_abundance = reference_prots[protID]
        if linear:
            # Use scaled protein abundance(0<abundance<1) for linear MOPA constraint to overcome the bias of simulation
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
            # option 2: Use the absolute protein abundance is better for quadratic MOPA
            ref_abundance = ref_abundance * prot.scaling_factor
            const = prob.Constraint(
                prot.concentration - dist,
                lb=ref_abundance,
                ub=ref_abundance,
                name="mopa_constraint_" + prot.id,
            )
            to_add.extend([dist, const])
            obj_vars.append(dist ** 2)
    model.add_cons_vars(to_add)

    # set MOPA objective
    if linear:
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
    else:
        model.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)


def add_mopa_ecGEM(model, reference_prots, linear: bool = True):
    """
    Add MOPA constraints and objective representing to the ecGEM model.
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
    for protID in reference_prots.index:
        prot=model.reactions.get_by_id(protID)
        ref_abundance = reference_prots[protID]
        if linear:
            components = sutil.add_absolute_expression(
                model,
                prot.flux_expression,
                name="mopa_dist_" + prot.id,
                difference=ref_abundance,
                add=False,
            )
            to_add.extend(components)
            obj_vars.append(components.variable)
        else:
            dist = prob.Variable("mopa_dist_" + prot.id)
            const = prob.Constraint(
                prot.flux_expression - dist,
                lb=ref_abundance,
                ub=ref_abundance,
                name="mopa_constraint_" + prot.id,
            )
            to_add.extend([dist, const])
            obj_vars.append(dist**2)
    model.add_cons_vars(to_add)
    if linear:
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
    else:
        model.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)
