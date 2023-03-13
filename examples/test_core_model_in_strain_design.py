# newly added
import sys
import os
sys.path.append('/Users/xluhon/Documents/GitHub/yetfl/code')
sys.path.append('/Users/xluhon/Documents/GitHub/etfl')
os.chdir('/Users/xluhon/Documents/GitHub/ETFLdesigner/examples')


from collections import namedtuple
import pandas as pd
import numpy  as np

from etfl.io.json import load_json_model, save_json_model
solver = 'optlang-gurobi'
ecoli = load_json_model("models/ecoli/ecoli_core_curated.json", solver=solver)
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

# flux variability analysis under the above max production
# before analysis, fix the bounds of production rate
ecoli.reactions.get_by_id('EX_succ_e').bounds = sol2.objective_value*0.8, sol2.objective_value

from pytfa.analysis import  variability_analysis
from etfl.optim.variables import mRNAVariable, EnzymeVariable
variables = EnzymeVariable
eva = variability_analysis(ecoli, variables)
eva_rxn = variability_analysis(ecoli, kind='reactions')


