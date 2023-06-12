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
    model.constraints["MODC_enzyme_fix_" + select_enzyme].ub = ub0
    model.constraints["MODC_enzyme_fix_" + select_enzyme].lb = ub0 # new
    return model



solver = 'optlang-gurobi'


# test related to strain design method development
ecoli = load_json_model("examples/models/ecoli/ecoli_core_curated.json", solver=solver)

# turn off by-product
ecoli.reactions.get_by_id('EX_for_e').bounds = 0,  0
ecoli.reactions.get_by_id('EX_ac_e').bounds = 0,  0
ecoli.reactions.get_by_id('EX_pyr_e').bounds = 0,  0
ecoli.reactions.get_by_id('EX_lac__D_e').bounds = 0,  0
ecoli.reactions.get_by_id('EX_etoh_e').bounds = 0,  0
ecoli.reactions.get_by_id('EX_acald_e').bounds = 0,  0

# definition of the initial state:
# firstly set growth as the objective function
ecoli.objective = "Biomass_Ecoli_core"
ecoli.objective_direction = 'max'
growth = ecoli.optimize()

# fix growth at 80% of its max value
ecoli.reactions.get_by_id('Biomass_Ecoli_core').bounds = growth.objective_value*0.8, growth.objective_value*0.8


# set the metabolite production（succinate） as the objective function
ecoli.objective = 'EX_succ_e'
ecoli.objective.direction = 'max'
sol2 = ecoli.optimize()
#ecoli.reactions.get_by_id('EX_succ_e').bounds = sol2.objective_value, sol2.objective_value
ecoli.reactions.get_by_id('EX_succ_e').bounds = 1, 1


# minimize glucose uptake
ecoli.objective = 'EX_glc__D_e'
ecoli.objective.direction = 'max'
glucose_uptake = ecoli.optimize()
# fix glucose uptake. Here actually the glucose uptake rate was fixed at a sub-optimal value
ecoli.reactions.get_by_id('EX_glc__D_e').bounds = glucose_uptake.objective_value, glucose_uptake.objective_value


# minimization the enzyme usage to get the protein abundance:
obj_expr = symbol_sum([ecoli.enzymes.dummy_enzyme.variable])
set_objective(ecoli, obj_expr)
ecoli.objective_direction = 'max'
ecoli.optimize()



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
out = prep_sol(substrate_uptake=-10, model=ecoli)
out_new = dict()
for x,y in out.items():
    if "EZ_" in x:
        x_new = x.replace('EZ_', '')
        #if "dummy_enzyme" in x_new:
        #    print(x_new)
        #else:
        out_new[x_new] = y
all_enzyme_ID = out_new.keys()




## give a disturbation
ecoli.reactions.get_by_id('Biomass_Ecoli_core').bounds = 0, 3 # first relax the above constraints
ecoli.reactions.get_by_id('EX_succ_e').bounds = 0, 15
# Here, we select PGI as a test example, first get its minimal abundance

# target enzyme

test_enz = 'ENO' # test another enzyme
test_enz = 'GLUDy' # test another enzyme
test_enz = 'PGI' # good example
test_enz = 'PGM_inferred_2' # good example
test_pro_abundance = out_new[test_enz]
ecoli = constrain_enzymes_based_abs_abundance(model=ecoli, select_enzyme='EZ_' + test_enz, ub0=test_pro_abundance*3)
## after disturbation, we need refix the max glucose uptake rate
ecoli.reactions.get_by_id('EX_glc__D_e').bounds = -5, -5






# run MOMA of protein resource adjust
def symbol_sum_MOMA(variables):
    from sympy import Add
    return Add(*variables)


"""
#  Not used
ss2 = []
for id in all_enzyme_ID:
    xx = (ecoli.enzymes.get_by_id(id).variable - out_new[id])**2
    ss2.append(xx)
"""

"""
# version 1
# when using this procedure, glucose uptake rate should be decreased to 4.25 mmol/gDW.h
ss2 = []
for id in all_enzyme_ID:
    xx = (1-(ecoli.enzymes.get_by_id(id).variable + 1e-07)/(out_new[id] + 1e-07))**2
    ss2.append(xx)
"""



# version 2 filter enzymes with zero values
out_new_filter = {}
for x, y in out_new.items():
    if y > 1e-07:
        out_new_filter[x] = y
all_enzyme_ID = out_new_filter.keys()
ss2 = []
for id in all_enzyme_ID:
    #xx = (ecoli.enzymes.get_by_id(id).variable - out_new[id])**2
    xx = (1-(ecoli.enzymes.get_by_id(id).variable)/(out_new[id]))**2
    ss2.append(xx)




"""
# this can't be used. ValueError: The given objective is invalid. Must be linear or quadratic.
ss = []
for id in all_enzyme_ID:
    xx = abs(ecoli.enzymes.get_by_id(id).variable - out_new[id])
    ss.append(xx)
"""

ecoli.objective = symbol_sum_MOMA(ss2)
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
result.to_excel("examples/result/abundance_before_and_after_MOMA.xlsx")


# plot
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


plt.figure()
sns.regplot(x=result['before'], y=result['after'], fit_reg=False)
sns.regplot(x=np.log10(result['before']), y=np.log10(result['after']), fit_reg=False)
plt.xlim(-10, 0)
plt.ylim(-10, 0)
plt.xlabel("log10(before)")
plt.ylabel("log10(after)")




