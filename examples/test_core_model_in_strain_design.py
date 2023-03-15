# newly added
import sys
import os
sys.path.append('/Users/xluhon/Documents/GitHub/yetfl/code')
sys.path.append('/Users/xluhon/Documents/GitHub/etfl')
os.chdir('/Users/xluhon/Documents/GitHub/ETFLdesigner/examples')


from etfl.io.json import load_json_model, save_json_model

from pytfa.optim.utils import symbol_sum
from cobra.util.solver import set_objective
solver = 'optlang-gurobi'
ecoli = load_json_model("models/ecoli/ecoli_core_curated.json", solver=solver)
# definition of the initial state:
# firstly set growth as the objective function
ecoli.objective = "Biomass_Ecoli_core"
ecoli.objective_direction = 'max'
growth = ecoli.optimize()

# evaluate the maximal growth under different glucose uptake rate
substrate = range(1, 15)
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
ecoli = load_json_model("models/ecoli/ecoli_core_curated.json", solver=solver)
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
















