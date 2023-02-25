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

#ecoli = load_json_model('../models/ME_ecoli_coreModel.json', solver=solver)
ecoli = load_json_model('../models/ecoli_core.json', solver=solver)



growth_reaction_id = 'Biomass_Ecoli_core'

glucose_uptake = -1
ecoli.reactions.EX_glc__D_e.lower_bound = glucose_uptake
ecoli.reactions.EX_glc__D_e.upper_bound = 0
# issue in growth coupling

# get GC_Biomass_Ecoli_core constraint
gc_growth=ecoli.constraints["GC_Biomass_Ecoli_core"]
ecoli.constraints["GC_Biomass_Ecoli_core"].ub= 3.0

# change upper bound of the constant allocation
pro_ratio = 0.051248
xx = ecoli.constant_allocation.enzyme_fix
xx.kwargs = {'ub': pro_ratio*0.75, 'lb': pro_ratio*0.75}
ecoli.constant_allocation.enzyme_fix = xx
cc = ecoli.constant_allocation.enzyme_fix
cc.constraint.ub = pro_ratio*0.75
cc.constraint.lb = pro_ratio*0.75

ecoli.constant_allocation.enzyme_fix = cc
"""
pp = ecoli.constant_allocation.prot_fix
pp.kwargs = {'ub': pro_ratio, 'lb': 0}
ecoli.constant_allocation.prot_fix = xx
mm = ecoli.constant_allocation.prot_fix
mm.constraint.lb = 0
mm.constraint.ub = pro_ratio
ecoli.constant_allocation.prot_fix = cc
"""
sol = ecoli.optimize()






# output fluxes and enzyme concentration
# how to output the pipetide concentration


def prep_sol(substrate_uptake, model, GLC_RXN_ID='r_1714'):
    from etfl.analysis.utils import enzymes_to_peptides_conc
    ret = {'obj':model.solution.objective_value,
           'mu':model.solution.fluxes.loc[model.growth_reaction.id],
           'available_substrate':-1*substrate_uptake,
           'uptake':-1*model.solution.fluxes[GLC_RXN_ID]
           }
#    for exch in model.exchanges:
#        ret[exch.id] = model.solution.fluxes.loc[exch.id]
    for rxn in model.reactions:
        ret[rxn.id] = model.solution.fluxes.loc[rxn.id]

    ret2 = dict()
    for enz in model.enzymes:
        ret2['EZ_' + enz.id] = model.solution.raw.loc['EZ_' + enz.id]
    peptides_conc = enzymes_to_peptides_conc(model, ret2)
    out = {**ret, **dict(peptides_conc)}
    return pd.Series(out)

solution = prep_sol(substrate_uptake=glucose_uptake, model=ecoli, GLC_RXN_ID='EX_glc__D_e')
solution2 = pd.DataFrame(solution)
solution2.columns =['value']
solution2.to_excel("../outputs/solution_for_core_ecoli.xlsx")









