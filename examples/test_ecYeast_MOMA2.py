import sys
import os
#os.chdir(r'D:\code\github\etfl\ETFLdesigner')
sys.path.append(r"/Users/xluhon/Documents/GitHub/ETFLdesigner/ETFLdesigner")


# load packages
from cobra.io import load_matlab_model, read_sbml_model
from pytfa.optim.utils import symbol_sum
from cobra.util.solver import set_objective
import pandas as pd
from ETFLdesigner.simulation.ecYeastFlux import *
from ETFLdesigner.manipulation.mainFunction import *



# heme model
ecYeast = read_sbml_model('examples/models/yeast/heme_ecYeastGEM.xml')
easy_model = produceRxnList(ecYeast)


product_name='heme_a'
product_id='EX_heme_a'

model2 = ecYeast.copy()
model2.reactions.get_by_id("r_1714_REV").bounds = (0, 100)  # D-glucose exchange (reversible)
model2.objective = {model2.reactions.r_2111: 1}
fba_solution1 = model2.optimize()
model2.objective = {model2.reactions.EX_heme_a: 1}
"""
fba_solution2 = model2.optimize()
result2 = productionEnvolpe(model=model2, biomassRxn='r_2111', targetRxn='EX_heme_a',substrateRxn='r_1714_REV', modelType = 'EC')
# plot
import matplotlib.pyplot as plt # For python 3.6: sudo apt-get install python3.6-tk   TO INSTALL Tkinter
x1 = result2['biomass']
y1l= result2['p_lb']
y1 = result2['p_ub']
plt.figure(figsize=(4,3))
plt.plot(x1,y1, 'k',  color='red')
plt.plot(x1,y1l, 'k',  color='red')
#plt.plot((0.32, 0.32), (0, 0.35518), 'k-', color='red')
plt.xlabel('growth(/h)')
plt.ylabel('production(mmol/gDW.h)')
plt.show()
"""

# test MOMA based on enzyme abundance adjustment
rxnID =[rxn.id for rxn in model2.reactions]
rxnID = [x for x in rxnID if "draw_prot_" in x]

# prepare the geneID
geneID = []
for x in ecYeast.genes:
    print(x.id)
    geneID.append(x.id)

# prepare ecModel for calculation
def RefSimulation(model0, D0, prxn, qp):
    """
    This funcion is used to simulate the chemostat growth of yeast
    Actually this function is general to solve the ecGEMs
    :param model0: a ecGEMs
    :param D0: a growth rate
    :return: solution of fluxes
    """
    growth = D0
    with model0:
        # model0 = ecYeastMinimalMedia(model0)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth, growth)
        model0.reactions.get_by_id(prxn).bounds = (qp, qp)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 1000)  # open the glucose
        model0.objective = {model0.reactions.r_1714_REV: -1}
        solution2 = model0.optimize()
        GR = solution2.fluxes["r_1714_REV"]  # get the glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (GR, GR * 1.001)
        model0.objective = {model0.reactions.prot_pool_exchange: -1}
        solution3 = model0.optimize()
    return solution3
# run MOMA of protein resource adjust
def symbol_sum_MOMA(variables):
    from sympy import Add
    return Add(*variables)

model2 = ecYeast.copy()
# set tolerance
model2.tolerance = 1e-9
model2.reactions.get_by_id("r_1714_REV").bounds = (0, 100)
fba_solution = RefSimulation(model0=model2, D0=0.3, prxn='EX_heme_a', qp=0.002)
qs = fba_solution['r_1714_REV']

# get the reference condition
uniprotID = 'P08417'
proID = 'draw_prot_' + uniprotID
target_pro = fba_solution[proID]




mutModel = model2.copy()
mutModel.reactions.get_by_id(proID).bounds = (target_pro*5, float('inf'))
test = mutModel.optimize()


# try to rewrite moma for ecModels specially!!
"""Provide minimization of metabolic adjustment (MOMA)."""
from typing import TYPE_CHECKING, Optional
from optlang.symbolics import Zero, add
from cobra.util import solver as sutil
from cobra.flux_analysis.parsimonious import pfba
from cobra.core import Model, Solution

mutModel.solver = sutil.choose_solver(mutModel, qp=True)
prob = mutModel.problem
v = prob.Variable("moma_old_objective")
# obj=test.objective_value
# v.bounds =obj,obj
c = prob.Constraint(
    mutModel.solver.objective.expression - v,
    lb=0.0,
    ub=0.0,
    name="moma_old_objective_constraint",
)
to_add = [v, c]
mutModel.objective = prob.Objective(Zero, direction="min", sloppy=True)
obj_vars = []
prot_rxns=mutModel.reactions.query(lambda x: 'draw_prot_' in x.id)
for r in prot_rxns:
    if 'draw_prot_' in r.id:
        flux = fba_solution.fluxes[r.id]
        print(flux)
        if flux > 0:
            print(r.id)
            components = sutil.add_absolute_expression(
                mutModel,
                r.flux_expression,
                name="moma_dist_" + r.id,
                difference=flux,
                add=False,
            )
            to_add.extend(components)
            obj_vars.append(components.variable)
mutModel.add_cons_vars(to_add)
mutModel.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})



# solve
fba_solution.fluxes['r_1714_REV']
fba_solution.fluxes['r_2111']
gluc=fba_solution.fluxes['r_1714_REV']
with mutModel:
    mutModel.reactions.get_by_id("r_1714_REV").bounds = (gluc, gluc)
    mutant = mutModel.optimize()
    mutant.fluxes['r_1714_REV']
    mutant.fluxes['r_2111']


# plot the result
mutant_df = mutant.fluxes
mutant_df = pd.DataFrame(mutant_df)
mutant_df = mutant_df[mutant_df.index.str.contains('draw_prot')]
mutant_df.columns =['mutant_flux']

wild_df = fba_solution.fluxes
wild_df = pd.DataFrame(wild_df)
wild_df = wild_df[wild_df.index.str.contains('draw_prot')]
wild_df.columns =['wild_flux']

merged = mutant_df.merge(wild_df, left_index=True, right_index=True, how='inner')



# plot
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

plt.figure()
sns.regplot(x=merged['wild_flux'], y=merged['mutant_flux'], fit_reg=False)
sns.regplot(x=np.log10(merged['wild_flux']), y=np.log10(merged['mutant_flux']), fit_reg=False)
plt.xlim(-10, 0)
plt.ylim(-10, 0)
plt.xlabel("log10(Measured_protein_level)")
plt.ylabel("log10(Predicted_protein_usage)")






# check the protein abundance change fold
import pandas as pd
prot_change=pd.Series()
for r in prot_rxns:
    if 'draw_prot_' in r.id:
        flux1 = fba_solution.fluxes[r.id]
        flux2= mutant.fluxes[r.id]
        if flux1 > 0:
            fc=flux2/flux1
            prot_change[r.id]=fc
        elif flux1 == 0:
            if flux2 > 0:
                fc=1000
            elif flux2 == 0:
                fc=1
            prot_change[r.id] = fc
        print(f"old_{r.id}",flux1,f'mut_{r.id}: ',flux2)
print('change more than 5 times:',len(prot_change[prot_change>5]))
print('change less than 0.2 times:',len(prot_change[prot_change<0.2]))
prot_change.to_excel('examples/result/test_moma.xlsx')


