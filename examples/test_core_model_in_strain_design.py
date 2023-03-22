# newly added
import sys
import os
sys.path.append('/Users/xluhon/Documents/GitHub/yetfl/code')
sys.path.append('/Users/xluhon/Documents/GitHub/etfl')
os.chdir('/Users/xluhon/Documents/GitHub/ETFLdesigner')


from etfl.io.json import load_json_model, save_json_model

from pytfa.optim.utils import symbol_sum
from cobra.util.solver import set_objective
import pandas as pd

# function
def prep_sol(substrate_uptake, model):
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
    return ret #pd.Series(ret)

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
                         lb=0,  # cannot be negative
                         ub=0) # once this value is changed, it can't be replaced by the new value
    ecoli.constraints["MODC_enzyme_fix_" + select_enzyme].ub = ub0
    ecoli.constraints["MODC_enzyme_fix_" + select_enzyme].lb = ub0 # new
    return model



solver = 'optlang-gurobi'


# test related to strain design method development
ecoli = load_json_model("examples/models/ecoli/ecoli_core_curated.json", solver=solver)



# definition of the initial state:
# firstly set growth as the objective function
ecoli.objective = "Biomass_Ecoli_core"
ecoli.objective_direction = 'max'
growth = ecoli.optimize()
# fix growth at 50% of its max value
ecoli.reactions.get_by_id('Biomass_Ecoli_core').bounds = growth.objective_value*0.5, growth.objective_value*0.5

# set the metabolite production（succinate） as the objective function
ecoli.objective = 'EX_succ_e'
ecoli.objective.direction = 'max'
sol2 = ecoli.optimize()

# minimization the enyzme usage to get the protein abundance:
ecoli.reactions.get_by_id('EX_succ_e').bounds = sol2.objective_value, sol2.objective_value
obj_expr = symbol_sum([ecoli.enzymes.dummy_enzyme.variable])
set_objective(ecoli, obj_expr)
ecoli.objective_direction = 'max'
ecoli.optimize()
sol3 = ecoli.optimize()

# how to output the enzyme concentration
ret = dict()
for enz in ecoli.enzymes:
    ret['EZ_' + enz.id] = ecoli.solution.raw.loc['EZ_' + enz.id]
ret['EZ_dummy_enzyme']
ez_rescale = dict()
for i, x in ret.items():
    print(i,x)
    i0 = i.replace('EZ_','')
    x_new = x*ecoli.enzymes.get_by_id(i0).scaling_factor
    ez_rescale[i] = x_new

pro_abundance = pd.Series(ez_rescale)
result = pd.DataFrame({'enzymeID': pro_abundance.index, 'reference': pro_abundance.values})
result.to_excel("examples/result/core_model_abundance.xlsx")

# obtain the reference protein abundance
GLC_RXN_ID ='EX_glc__D_e'

out = prep_sol(substrate_uptake=-10, model=ecoli)
out_new = dict()
for x,y in out.items():
    if "EZ_" in x:
        x_new = x.replace('EZ_', '')
        out_new[x_new] = y
all_enzyme_ID = out_new.keys()







# give a disturbation
ecoli.reactions.get_by_id('Biomass_Ecoli_core').bounds = 0, growth.objective_value*0.5# first relax the above constraints
ecoli.reactions.get_by_id('EX_succ_e').bounds = sol2.objective_value, 1000
#ecoli = constrain_enzymes_based_abs_abundance(model=ecoli, select_enzyme = 'EZ_TPI', ub0=0.005)
ecoli = constrain_enzymes_based_abs_abundance(model=ecoli, select_enzyme = 'EZ_PGI', ub0=0.00173*5)


# run MOMA of protein resource adjust
def symbol_sum_MOMA(variables):
    from sympy import Add
    return Add(*variables)
ss = []
for id in all_enzyme_ID:
    xx = (ecoli.enzymes.get_by_id(id).variable - out_new[id])**2
    ss.append(xx)

"""
# this can't be used. ValueError: The given objective is invalid. Must be linear or quadratic.
ss = []
for id in all_enzyme_ID:
    xx = abs(ecoli.enzymes.get_by_id(id).variable - out_new[id])
    ss.append(xx)
"""
ecoli.objective = symbol_sum_MOMA(ss)
ecoli.objective_direction = 'min'
sol4 = ecoli.optimize()
out2 = prep_sol(substrate_uptake=-10, model=ecoli)
out_new2 = dict()
for x,y in out2.items():
    if "EZ_" in x:
        x_new = x.replace('EZ_', '')
        out_new2[x_new] = y

df = pd.Series(out_new)
df2 = pd.Series(out_new2)

result = pd.concat([df, df2], axis=1)
result.columns =['before','after']



# plot
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

plt.figure()
sns.regplot(x=result['before'], y=result['after'], fit_reg=False)
sns.regplot(x=np.log10(result['before']), y=np.log10(result['after']), fit_reg=False)
plt.xlim(-10, -3)
plt.ylim(-10, -3)
plt.xlabel("log10(before)")
plt.ylabel("log10(after)")






