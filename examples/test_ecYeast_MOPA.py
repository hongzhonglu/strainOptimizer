# -*- coding: utf-8 -*-
# date : 2024/3/12 
# author : wangh
from optlang.symbolics import Zero, add
from cobra.util import solver as sutil
from strainOptimizer.simulation.pprotFBA import pprotFBA_prot_conc
import pandas as pd
import os


def mopa(model, reference_prots, linear: bool = True) :
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
        

def compare_mutant_vs_reference(ref,mutant,threshold=5):
    """
    Compare the protein abundace change of mutant vs reference
    """
    prot_change = pd.Series()
    for r in ref.index:
        if 'draw_prot_' in r:
            flux1 = ref[r]
            flux2 = mutant.fluxes[r]
            if flux1 > 0:
                fc = flux2 / flux1
                prot_change[r] = fc
            elif flux1 == 0:
                if flux2 > 0:
                    fc = 1000
                elif flux2 == 0:
                    fc = 1
                prot_change[r] = fc
    print('change more than 5 times:', len(prot_change[prot_change > threshold]))
    print('change less than 0.2 times:', len(prot_change[prot_change < 1/threshold]))
    return prot_change



if __name__ == "__main__":
    from strainOptimizer.io import load_model

    # set working directory
    os.chdir(r'D:\code\github\strainOptimizer')

    # load model
    ecYeast = load_model("examples/models/yeast/heme_ecYeastGEM.xml",model_type='ecGEM')

    product_name = 'heme_a'
    product_id = 'EX_heme_a'
    growth='r_2111'
    c_source='r_1714_REV'
    c_uptake=10

    # set wild type condition as the reference
    ref_prots = pprotFBA_prot_conc(model=ecYeast,
                                   targetID=growth,
                                   c_source=c_source,
                                   c_uptake=c_uptake,
                                   model_type='ecGEM')


    # select overexpression protein as test example
    uniprotID = 'P08417'
    proID = 'draw_prot_' + uniprotID

    ref_prot_conc = ref_prots[proID]

    mutModel = ecYeast.copy()
    mutModel.reactions.get_by_id(proID).bounds = (ref_prot_conc * 5, float('inf'))

    # test MOPA
    sol = mopa(model=mutModel, reference_prots=ref_prots, linear=False)

    sol2= mopa(model=mutModel, reference_prots=ref_prots, linear=True)


    # compare the protein abundance change fold
    mopa_change = compare_mutant_vs_reference(ref=ref_prots,mutant=sol,threshold=5)   # quadratic optimization is more suitable for ecGEM

    mopa_linear_change = compare_mutant_vs_reference(ref=ref_prots,mutant=sol2,threshold=2)   # linear MOPA shows too less proteins with changed abundance.


