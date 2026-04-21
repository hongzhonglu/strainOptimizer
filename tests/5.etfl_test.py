# -*- coding: utf-8 -*-
from strainOptimizer.io import load_model
import numpy as np
import cobra
from strainOptimizer.simulation.pprotFBA import ppFBA
# cobra.Configuration().tolerance=1e-9

model=load_model(filename='examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json',
                 solver='optlang-gurobi',
                 model_type='etfl')

c_uptake=1
c_source="r_1714"      # glucose exchange rxn
target_id='r_1589'  # 2-phenylethanol exchange rxn

model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake

model.objetive=model.growth_reaction
model.slim_optimize()

growth_values=np.linspace(0.02,0.08,10)
for growth_value in growth_values:
    # growth_value=round(growth_value,4)
    with model:
        model.reactions.get_by_id(model.growth_reaction.id).bounds=growth_value*0.99,growth_value
        model.objective=target_id
        # solution=model.optimize()
        solution=ppFBA(model,target_id,c_source,c_uptake=c_uptake,model_type='etfl')
        print('growth_value: %s, solution: %s'%(growth_value,solution))
