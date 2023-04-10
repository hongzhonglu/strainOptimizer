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


# how to add rxns
# with Itaconic acid as an example
from cobra import  Reaction, Metabolite
yeast.add_metabolites(Metabolite('s_4349_c', compartment='c', formula='C5H4O4', name='itaconate'))
yeast.add_metabolites(Metabolite('s_4350_e', compartment='e', formula='C5H4O4', name='itaconate'))

"""
# check the existence of metabolites
for met in yeast.metabolites:
    print(met.id)
    if met.id.find('s_0794') != -1:
        print(met.id)"""

# define the two reactions
dict2 = {yeast.metabolites.get_by_id('s_0516_c'): -1,
         #yeast.metabolites.get_by_id('s_0794_c'): -1, # no H+ in the yetfl model, why?
         yeast.metabolites.get_by_id('s_4349_c'): 1,
         yeast.metabolites.get_by_id('s_0456_c'): 1,
        }


# reaction itaconate[c] => itaconate[e]
dict3 = {yeast.metabolites.get_by_id('s_4349_c'): -1,
         yeast.metabolites.get_by_id('s_4350_e'): 1,
        }

#  sink reaction
dict4 = {yeast.metabolites.get_by_id('s_4350_e'): -1
        }


# add the reactions into the model
# 2
# using the general procedures to add the new function
reaction = Reaction('new_cis_aconitate_decarboxylase')
reaction.name = 'new_cis_aconitate_decarboxylase'
reaction.subsystem = 'new added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict2)
reaction.gene_reaction_rule = 'AtCAD'
yeast.add_reactions([reaction])


# 3
reaction = Reaction('new_itaconate_trasnport')
reaction.name = 'new_itaconate_trasnport'
reaction.subsystem = 'new added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict3)
reaction.gene_reaction_rule = ''
yeast.add_reactions([reaction])

# 4
reaction = Reaction('exchange_itaconate')
reaction.name = 'exchange_itaconate'
reaction.subsystem = 'new added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict4)
reaction.gene_reaction_rule = ''
yeast.add_reactions([reaction])



# solve
yeast.objective = 'exchange_itaconate'
yeast.objective.direction = 'max'
yeast.reactions.r_4041.lower_bound = 0.2
sol = yeast.optimize()
sol.fluxes['r_1714']

print(yeast.summary())


