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

def constrain_enzymes_based_abs_abundance(model, select_enzyme=['EZ_TPI'], ub0=1):
    from etfl.optim.variables import EnzymeVariable
    from etfl.optim.constraints import ModelConstraint
    # a function to add a constraint on total amount of enzymes based on their
    # fraction from total amount of proteins (should be before adding dummy)
    enz_vars = model.get_variables_of_type(EnzymeVariable)
    #enz_vars = [x for x in enz_vars if x.name not in exclusion]
    for en in select_enzyme:
        print(en)
        single_vars = [x for x in enz_vars if x.name == en]
        expr = symbol_sum([x for x in single_vars])
        model.add_constraint(kind=ModelConstraint,
                             hook=model,
                             expr=expr,
                             id_='enzyme_fix_' + en,
                             lb=0,  # cannot be negative
                             ub=0)  # once this value is changed, it can't be replaced by the new value
        model.constraints["MODC_enzyme_fix_" + en].ub = ub0
        model.constraints["MODC_enzyme_fix_" + en].lb = ub0  # new

    return model

solver = 'optlang-gurobi'


# test related to strain design method development
ecoli = load_json_model("examples/models/ecoli/ecoli_core_curated.json", solver=solver)



# definition of the initial state:
# firstly set growth as the objective function

# turn off by-product
def solveModel(model=ecoli):
    import math
    model.reactions.get_by_id('EX_for_e').bounds = 0, 0
    model.reactions.get_by_id('EX_ac_e').bounds = 0, 0
    model.objective = "Biomass_Ecoli_core"
    model.objective_direction = 'max'
    #growth = model.optimize()
    growth = model.slim_optimize()
    if math.isnan(growth):
        qpmax = None
        qpmin = None
    else:
        # fix growth at 50% of its max value
        model.reactions.get_by_id('Biomass_Ecoli_core').bounds = 0.2, 3
        model.reactions.get_by_id('EX_succ_e').bounds = 0, 100
        # fix growth at its max
        model.reactions.get_by_id(
            'Biomass_Ecoli_core').bounds = growth * 0.95, growth * 0.95
        # set the metabolite production（succinate） as the objective function. Then minimize the production:
        model.objective = 'EX_succ_e'
        model.objective.direction = 'max'
        qpmax = model.slim_optimize()
        model.objective.direction = 'min'
        qpmin = model.slim_optimize()
    return growth, qpmax, qpmin








# minimization the enzyme usage to get the protein abundance:
ecoli.reactions.get_by_id('EX_succ_e').bounds =2.74*0.5, 2.74*0.5
obj_expr = symbol_sum([ecoli.enzymes.dummy_enzyme.variable])
set_objective(ecoli, obj_expr)
ecoli.objective_direction = 'max'
ecoli.optimize()
sol3 = ecoli.optimize()


# get enzyme list which need to be tuned
out2 = prep_sol(substrate_uptake=-10, model=ecoli)
out_new2 = dict()
for x,y in out2.items():
    if "EZ_" in x:
        x_new = x.replace('EZ_', '')
        out_new2[x_new] = y
df2 = pd.Series(out_new2)
enzyme = list(df2.index)
list1 = df2.to_list()
list_ref = [1] * len(list1)
index_ref = list(range(0,len(list1)))

import numpy as np
np.random.choice(index_ref, 3, replace=False)
total_sample = 1000
result = []
for i in range(0, total_sample):
    print(i)
    s = list(np.random.choice(index_ref, 3, replace=False))
    result.append(s)

first_generation = []
for item in result:
    print(item)
    exist = []
    for ii, xx in enumerate(list_ref):
        if ii in item:
            exist.append(0)
        else:
            exist.append(xx)
    first_generation.append(exist)

# calculate fitness of strain
fitness = []
for ii, xx in enumerate(first_generation):
    print(ii, xx)
    ss = xx
    delete_gene = [x for x,y in zip(enzyme,xx) if y==0]
    delete_gene = ['EZ_'+ x for x in delete_gene]
    # update the model
    new_model = constrain_enzymes_based_abs_abundance(model=ecoli.copy(), select_enzyme=delete_gene, ub0=0)
    # solve the new model
    growth, max, min = solveModel(model=new_model)
    out =[growth, max, min]
    fitness.append(out)
    print(growth,max, min)













# define the max production rate
# fix growth at its max
ecoli.reactions.get_by_id('Biomass_Ecoli_core').bounds = 0.2, 0.2
ecoli.reactions.get_by_id('EX_succ_e').bounds = 0, 1000
# set the metabolite production（succinate） as the objective function. Then minimize the production:
ecoli.objective = 'EX_succ_e'
ecoli.objective.direction = 'max'
qpmax = ecoli.optimize()





















