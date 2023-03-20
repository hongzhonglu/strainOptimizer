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


solver = 'optlang-gurobi'
ecoli = load_json_model("examples/models/ecoli/ecoli_core_curated.json", solver=solver)
# definition of the initial state:
# firstly set growth as the objective function
ecoli.objective = "Biomass_Ecoli_core"
ecoli.objective_direction = 'max'
growth = ecoli.optimize()

# evaluate the maximal growth under different glucose uptake rate
substrate = range(1, 10)
growth = []
for s in substrate:
    print(s)
    ecoli.reactions.get_by_id('EX_glc__D_e').bounds = -s, 0
    ecoli.optimize()
    sol3 = ecoli.optimize()
    v = sol3.objective_value
    growth.append(v)


# evaluate dummy enzyme usage under different growth
growth = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
dummy_pro_usage = []
for g in growth:
    print(g)
    ecoli.reactions.get_by_id('Biomass_Ecoli_core').bounds = g, g
    obj_expr = symbol_sum([ecoli.enzymes.dummy_enzyme.variable])
    set_objective(ecoli, obj_expr)
    ecoli.objective_direction = 'max'
    ecoli.optimize()
    sol3 = ecoli.optimize()
    dummy_pro = sol3.objective_value
    dummy_pro_usage.append(dummy_pro)





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



# conduct the flux variability analysis
# in the current analysis, it seems that it not consider the scaling factor:
from pytfa.analysis import  variability_analysis
from etfl.optim.variables import mRNAVariable, EnzymeVariable
variables = EnzymeVariable
eva = variability_analysis(ecoli, variables)
# merge it with the reference protein abundance
eva.to_excel("examples/result/enzyme_fva_core_model.xlsx")



















