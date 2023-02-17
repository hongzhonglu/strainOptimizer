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

uptake_range = pd.Series(list(range(10)))

def simulate(uptake, model):
    model.reactions.r_1714.lower_bound = -1*uptake
    model.reactions.r_1714.upper_bound = 0
    model.objective = growth_reaction_id
    model.objective.direction = 'max'
    try:
        sol = model.optimize()
    except Exception:
        return pd.Series([np.nan, np.nan])
    return pd.Series([sol.objective_value, model.reactions.r_1714.flux*-1])
"""
try:
    me_data = pd.read_csv('outputs/benchmark_T1E1N1.csv')
except FileNotFoundError:
    me_data = uptake_range.apply(simulate, args=[yeast])
    me_data.columns = ['mu','uptake']"""

me_data = []
for x in uptake_range:
    print(x)
    ss = simulate(x, yeast)
    me_data.append(ss)
