# newly added
import sys
import os
sys.path.append('/Users/xluhon/Documents/GitHub/yetfl/code')
sys.path.append('/Users/xluhon/Documents/GitHub/etfl')

# check and set the current directory
#os.chdir('/Users/xluhon/Documents/GitHub/ETFLdesigner')


from collections import namedtuple
import pandas as pd
import numpy  as np

from etfl.io.json import load_json_model, save_json_model

solver = 'optlang-gurobi'

# ec_cobra.reactions.ATPM.lower_bound = 0

#ecoli = load_json_model('../models/ME_ecoli_coreModel.json', solver=solver)
ecoli = load_json_model('examples/models/ecoli/ecoli_core.json', solver=solver)



growth_reaction_id = 'Biomass_Ecoli_core'

glucose_uptake = -10
ecoli.reactions.EX_glc__D_e.lower_bound = glucose_uptake
ecoli.reactions.EX_glc__D_e.upper_bound = 0
# issue in growth coupling

# get GC_Biomass_Ecoli_core constraint
gc_growth=ecoli.constraints["GC_Biomass_Ecoli_core"]
ecoli.constraints["GC_Biomass_Ecoli_core"].ub= 3.0

# change upper bound of the constant allocation
pro_ratio = 0.051248*1.5
xx = ecoli.constant_allocation.enzyme_fix
xx.kwargs = {'ub': pro_ratio, 'lb': pro_ratio}
ecoli.constant_allocation.enzyme_fix = xx
cc = ecoli.constant_allocation.enzyme_fix
cc.constraint.ub = pro_ratio
#cc.constraint.lb = pro_ratio # as total protein resource is fixed, here we not fix the lower bound
cc.constraint.lb = 0

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


# remove unused genes and then further test model function
gene_list = []
for xx in ecoli.transcription_reactions:
    print(xx.id)
    gene = xx.id.replace("_transcription", "")
    gene_list.append(gene)

#gene_list.append('b3739')
remove_gene = []
for gene in ecoli.genes:
    if gene.id not in gene_list:
        print(gene.id)
        remove_gene.append(gene.id)
for xx in remove_gene:
    ecoli.genes.remove(xx)

# save and test the model
save_json_model(ecoli, "examples/models/ecoli_core_curated.json")
ecoli = load_json_model("examples/models/ecoli_core_curated.json", solver=solver)
sol_test = ecoli.optimize()







# how to add constraint
from etfl.optim.constraints import ModelConstraint
from pytfa.optim.utils import symbol_sum


ecoli.objective = symbol_sum([ecoli.reactions.get_by_id(x.id).forward_variable + ecoli.reactions.get_by_id(x.id).reverse_variable for x in ecoli.reactions if x.id == 'EX_glc__D_e'])
ecoli.objective_direction = 'min'
ecoli.slim_optimize()


upt = 2
tol = 0.02
#expr0 = ecoli.objective.expression
expr0 = symbol_sum([ecoli.reactions.get_by_id(x.id).forward_variable + ecoli.reactions.get_by_id(x.id).reverse_variable for x in ecoli.reactions if x.id == 'EX_glc__D_e'])
sub_cons = ecoli.add_constraint(kind=ModelConstraint,
                                hook=ecoli,
                                expr=expr0,
                                id_='fix_substrate',
                                lb=upt - abs(tol * upt),
                                ub=upt + abs(tol * upt),)
# note: for this kind of constraint, MODC_fix_substrate, will be the name of constraint in the model
# check the constraints
# new_cons = ecoli.constraints['fix_substrate']
ecoli.objective = "Biomass_Ecoli_core"
ecoli.objective_direction = 'max'
sol = ecoli.optimize()
# if removing the constraints
#ecoli.remove_constraint(sub_cons)
#sol = ecoli.optimize()




# how to add enzyme absolute concentration as constraits
from etfl.optim.variables import EnzymeVariable
from etfl.optim.constraints import ConstantAllocation
def constrain_enzymes_based_abs_abundance(model=ecoli, select_enzyme='EZ_TPI', ub0=1):
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
    return model

# test the function
ecoli = constrain_enzymes_based_abs_abundance(model=ecoli, select_enzyme = 'EZ_TPI', ub0=0)
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
solution2.to_excel("examples/result/solution_for_core_ecoli.xlsx")


# predict the maximum yield of product
#  add a sink reaction
from cobra import Reaction
dict4 = {ecoli.metabolites.get_by_id('mal__L_c'): -1}
reaction = Reaction('exchange_malic_acid')
reaction.name = 'exchange_malic_acid'
reaction.subsystem = 'new added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict4)
reaction.gene_reaction_rule = ''
ecoli.add_reactions([reaction])
# firstly set growth as the objective function
ecoli.objective = "Biomass_Ecoli_core"
ecoli.objective_direction = 'max'
growth = ecoli.optimize()
# fix growth at 80% of its max value
ecoli.reactions.get_by_id('Biomass_Ecoli_core').bounds = growth.objective_value*0.8, growth.objective_value*0.8

# set the metabolite production（malic acid） as the objective function
ecoli.objective = 'exchange_malic_acid'
ecoli.objective.direction = 'max'
sol1 = ecoli.optimize()

# set the metabolite production（succinate） as the objective function
ecoli.objective = 'EX_succ_e'
ecoli.objective.direction = 'max'
sol2 = ecoli.optimize()


# conduct the flux variability analysis
# in the current analysis, it seems that it not consider the scaling factor:
from pytfa.analysis import  variability_analysis
from etfl.optim.variables import mRNAVariable, EnzymeVariable
from cobra.util.solver import set_objective
variables = EnzymeVariable
eva = variability_analysis(ecoli, variables)
# merge it with the reference protein abundance
eva.to_excel("examples/result/enzyme_fva_core_model.xlsx")


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
