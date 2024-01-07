'''add new enzyme into prot_pool constraint'''
from strainOptimizer.io import load_model
from cobra import Metabolite,Reaction

# load ecModel
model= load_model('examples/models/yeast/ecYeastGEM_batch.xml',model_type='ecGEM')

met_prot_pool=model.metabolites.get_by_id('prot_pool[c]')
exchange_prot_pool=model.reactions.get_by_id('prot_pool_exchange')

# 1. add new enzyme mets
new_enzyme=Metabolite('new_enzyme[c]',name='new enzyme',compartment='c')
model.add_metabolites(new_enzyme)
new_MW=70 # 1 kDa= 1000 g/mol = 1 g/mmol

# 2. add draw protein reaction to add new enzyme into prot_pool
new_draw_rxn=Reaction('draw_new_enzyme')
new_draw_rxn.add_metabolites({new_enzyme:1,
                              met_prot_pool:-new_MW})

model.add_reactions([new_draw_rxn])

