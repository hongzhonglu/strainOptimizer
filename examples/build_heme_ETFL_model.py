'''Build heme a production model in ETFL by add heme a exchange reaction'''
from strainOptimizer.io import load_model

# load model
model=load_model('examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json', solver='optlang-gurobi',model_type='etfl')

heme_metID='s_0811_c'
heme_met=model.metabolites.get_by_id(heme_metID)

# add heme a demand reaction to simulate heme a production
from cobra import Reaction
reaction=Reaction('EX_heme_a')
reaction.name='heme a production'
reaction.subsystem='heme a production'
reaction.lower_bound=0
reaction.upper_bound=1000
reaction.add_metabolites({heme_met:-1})
model.add_reactions([reaction])

# set the objective function
model.objective='EX_heme_a'
model.objective_direction='max'

# set the glucose uptake rate
c_source='r_1714'
c_uptake=1
model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake

# optimize the model
sol=model.optimize()

# save the model
from etfl.io.json import save_json_model
save_json_model(model,'examples/models/yeast/heme_cEFL.json')
