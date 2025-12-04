# newly added
import sys
import os
# sys.path.append('/Users/xluhon/Documents/GitHub/yetfl/code')
# sys.path.append('/Users/xluhon/Documents/GitHub/etfl')
# os.chdir('/Users/xluhon/Documents/GitHub/yetfl/code')





from collections import namedtuple
import pandas as pd
import numpy  as np

from etfl.io.json import load_json_model

solver = 'optlang-gurobi'
solver= 'optlang-cplex'


# ec_cobra.reactions.ATPM.lower_bound = 0
growth_reaction_id = 'r_4041'
yeast = load_json_model('examples/models/yeast/yeast8_cEFL_2584_enz_128_bins__20221031_130538.json', solver=solver)
# change solver
yeast.solver = solver
yeast.reactions.r_1714.lower_bound = -10
yeast.reactions.r_1714.upper_bound = 0
yeast.objective = growth_reaction_id
yeast.objective.direction = 'max'
sol = yeast.optimize()



# print rxn the model
for rxn in yeast.reactions:
    print(rxn.id)


out = yeast.solution.raw.to_frame()
#out.to_excel("../outputs/whole_output_cEFL.xlsx")








#############################################################
# how to output concentration
#############################################################

# here not sure whether need to rescale the output
# be careful: this function may not output the real concentration of mRNA OR enzyme!
def prep_sol(substrate_uptake, model, GLC_RXN_ID='r_1714'):

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
        ret['EZ_'+ enz.id] = model.solution.raw.loc['EZ_'+ enz.id]
    for mrna in model.genes:
        print(mrna.id)
        ret['MR_'+ mrna.id] = model.solution.raw.loc['MR_'+ mrna.id]
    for ribosome in model.peptides:
        print(ribosome.id)
        if ribosome.id == 'dummy_peptide':
            ribosome.id='dummy_gene'
            ret['RP_' + ribosome.id] = model.solution.raw.loc['RP_' + ribosome.id]
        else:
            ret['RP_'+ ribosome.id] = model.solution.raw.loc['RP_'+ ribosome.id]
    return pd.Series(ret)

all_solution = prep_sol(substrate_uptake=-10, model=yeast)



# how to output the pipetide concentration
from etfl.analysis.utils import enzymes_to_peptides_conc
ret = dict()
for enz in yeast.enzymes:
    ret['EZ_' + enz.id] = yeast.solution.raw.loc['EZ_' + enz.id]
peptides_conc = enzymes_to_peptides_conc(yeast, ret)

# how to output the enzyme concentration
ez_rescale = dict()
for i, x in ret.items():
    print(i,x)
    i0 = i.replace('EZ_','')
    x_new = x*yeast.enzymes.get_by_id(i0).scaling_factor
    ez_rescale[i] = x_new



#############################################################
# how to adjust the kcat of specific rxns?
#############################################################
yeast.reactions.r_1714.lower_bound = -2
yeast.reactions.r_1714.upper_bound = 0
yeast.objective = growth_reaction_id
yeast.objective.direction = 'max'
sol = yeast.optimize()

yeast.enzymes.get_by_id('YAL038W_enzyme').kcat_fwd = 100
yeast.enzymes.get_by_id('YAL038W_enzyme').kcat_bwd = 100

yeast.enzymes.get_by_id('YOR347C_enzyme').kcat_fwd = 100
yeast.enzymes.get_by_id('YOR347C_enzyme').kcat_bwd = 100

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


#############################################################
# how to print the result
#############################################################
from etfl.optim.utils import fix_growth, release_growth, \
                            get_active_growth_bounds, safe_optim
from etfl.optim.config import standard_solver_config
standard_solver_config(yeast, verbose=False)
yeast.optimize()

print('Objective            : {}'.format(yeast.solution.objective_value))
print(' - Glucose uptake    : {}'.format(yeast.reactions.r_1714.flux))
print(' - Growth            : {}'.format(yeast.growth_reaction.flux))
print(' - Ribosomes produced: {}'.format(yeast.ribosome))
print(' - RNAP produced: {}'.format(yeast.rnap))
sol_max = safe_optim(yeast)









#############################################################
# flux variability analysis
#############################################################
from etfl.io.json import load_json_model

from etfl.optim.variables import mRNAVariable, EnzymeVariable
from etfl.optim.utils import fix_integers

from pytfa.analysis import  variability_analysis,           \
                            apply_reaction_variability,     \
                            apply_generic_variability,      \
                            apply_directionality

from etfl.analysis.utils import enzymes_to_peptides_conc

import pandas as pd



def print_sol(model):
    print('Objective            : {}'.format(model.solution.objective_value))
    #print(' - Glucose uptake    : {}'.format(model.solution.raw.loc['EX_glc__D_e']))
    #print(' - Growth            : {}'.format(model.solution.raw.loc[model.growth_reaction.id]))
    #print(' - Ribosomes produced: {}'.format(model.solution.raw.loc['EZ_rib']))
    #print(' - RNAP produced: {}'.format(model.solution.raw.loc['EZ_rnap']))
    try:
        print(' - DNA produced: {}'.format(model.solution.raw.DN_DNA))
    except AttributeError:
        pass



uptake_high = -12.5 #mmol/gDW/h

uptake_low  = -1 # mmol/gDW/h

# ecoli.growth_reaction.lower_bound = ecoli.solution.f - ecoli.solver.configuration.tolerances.feasibility
yeast.optimize()


yeast.solver.configuration.verbosity = 1
yeast.solver.configuration.tolerances.feasibility = 1e-9
try:
    yeast.solver.problem.Params.NumericFocus = 3
except AttributeError:
    pass

yeast.solver.configuration.presolve = True

# print('Fixing integers..')
# continuous_model = fix_integers(yeast)
continuous_model = yeast

# continuous_model.remove_constraint(continuous_model.sos1_constraint.interpolation_integer_SOS1)
print('optimizing..')
continuous_model.optimize()


def get_mu_bin(model, mu):
    closeness = [abs(x[0]-mu) for x in yeast.mu_bins]
    ix = closeness.index(min(closeness))
    return model.mu_bins[ix][1]

# Calculate variability analysis on all continuous variables

for glc_uptake in [uptake_high,uptake_low]:
    continuous_model.growth_reaction.lower_bound = 0
    continuous_model.growth_reaction.upper_bound = 10
    continuous_model.reactions.r_1714.lower_bound = glc_uptake - 0.01
    continuous_model.reactions.r_1714.upper_bound = glc_uptake + 0.01

    continuous_model.optimize()
    mu = continuous_model.growth_reaction.flux
    mu_lo, mu_hi = get_mu_bin(continuous_model, mu)

    continuous_model.growth_reaction.upper_bound = mu_hi
    continuous_model.growth_reaction.lower_bound = mu_lo

    print_sol(continuous_model)

    variables = EnzymeVariable

    eva = variability_analysis(continuous_model, variables)
    peptides_conc_min = pd.Series(enzymes_to_peptides_conc(continuous_model, eva['minimum']))
    peptides_conc_max = pd.Series(enzymes_to_peptides_conc(continuous_model, eva['maximum']))
    peptides_conc = pd.concat([peptides_conc_min,peptides_conc_max], axis=1)
    peptides_conc.to_csv('../outputs/iJO_vETFL_low_hi_{}_pep.csv'.format(glc_uptake))

    rescale = lambda row: yeast.enzymes.get_by_id(row.name[3:]).scaling_factor * row
    eva_real = eva.apply(rescale, axis=1)
    eva_real.to_csv('../outputs/iJO_vETFL_1783_low_hi_{}_enz.csv'.format(glc_uptake))

    # mva = variability_analysis(continuous_model, mRNAVariable)
    # rescale = lambda row: yeast.mrnas.get_by_id(row.name[3:]).scaling_factor * row
    # mva_real = mva.apply(rescale, axis=1)
    # mva_real.to_csv('outputs/iJO_vETFL_low_hi_{}_mrna.csv'.format(glc_uptake))



#############################################################
# test the time for the optimization
#############################################################
model.solver.configuration.timeout = 18000
model.solver.configuration.tolerances.optimality = 1e-4
model.solver.problem.Params.MIPFocus = 0
# growth_uptake_config(model)
model.warm_start = None







