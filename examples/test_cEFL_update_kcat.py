# newly added
import sys
import os
sys.path.append('/Users/xluhon/Documents/GitHub/yetfl/code')
sys.path.append('/Users/xluhon/Documents/GitHub/etfl')
os.chdir('/Users/xluhon/Documents/GitHub/yetfl/code')





from collections import namedtuple
import pandas as pd
import numpy  as np

from etfl.io.json import load_json_model

solver = 'optlang-gurobi'

# ec_cobra.reactions.ATPM.lower_bound = 0
growth_reaction_id = 'r_4041'
yeast = load_json_model('../models/yeast8_cEFL_2584_enz_128_bins__20221031_130538.json', solver=solver)
yeast.reactions.r_1714.lower_bound = -10
yeast.reactions.r_1714.upper_bound = 0
yeast.objective = growth_reaction_id
yeast.objective.direction = 'max'
sol = yeast.optimize()


#############################################################
# how to adjust the kcat of specific rxns?
#############################################################
yeast.reactions.r_1714.lower_bound = -2
yeast.reactions.r_1714.upper_bound = 0
yeast.objective = growth_reaction_id
yeast.objective.direction = 'max'
sol = yeast.optimize()

yeast.enzymes.get_by_id('YAL038W_enzyme').kcat_fwd = 10
yeast.enzymes.get_by_id('YAL038W_enzyme').kcat_bwd = 10

yeast.enzymes.get_by_id('YOR347C_enzyme').kcat_fwd = 10
yeast.enzymes.get_by_id('YOR347C_enzyme').kcat_bwd = 10

# update variable and constraints attributes?
rid = 'r_0962'
r = yeast.reactions.get_by_id(rid)
yeast.apply_enzyme_catalytic_constraint(r)
yeast._push_queue()
yeast.regenerate_constraints()
yeast.regenerate_variables()
sol2 = yeast.optimize()

constraints = yeast.constraints.FC_r_0962
yeast.remove_constraint(cons=constraints)
#yeast.regenerate_constraints() # no use to add back the constraint of r_0962
#yeast.repair() # no use to add back the constraint of r_0962
# thus we need other strategies
yeast.coupling_dict.update()
r = yeast.reactions.get_by_id('r_0962')
r.add_enzymes([yeast.enzymes.get_by_id('YOR347C_enzyme'), yeast.enzymes.get_by_id('YAL038W_enzyme')])
yeast.apply_enzyme_catalytic_constraint(r)
yeast._push_queue()
yeast.regenerate_variables()
yeast.regenerate_constraints()

constraints = yeast.constraints.FC_r_0962
print(constraints.expression)
sol3 = yeast.optimize()

