# -*- coding: utf-8 -*-
from strainOptimizer.io import load_model
import os
import pandas as pd
import numpy as np
import seaborn as sns

os.chdir(r'D:\code\github\strainOptimizer')
solver = 'optlang-gurobi'

# load model
modelParams_dict={
    'ecGEM':{'path':'examples/models/yeast/ecYeastGEM_batch.xml',
             'growth_id':'r_2111',
             'c_source':'r_1714_REV'
             },
    'etfl':{'path':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json',
           'growth_id':'r_2111',
           'c_source':'r_1714'
           },
    'gem':{'path':'examples/models/yeast/yeast-GEM.xml',
           'growth_id':'r_2111',
           'c_source':'r_1714'}
}

c_uptakes=np.linspace(0.5,15.5,25)
growth_sim={}
for type,modelParams in modelParams_dict.items():
    type='ecGEM'
    modelParams=modelParams_dict[type]
    model=load_model(modelParams['path'],model_type=type,solver=solver)
    growth_list=list()
    for c_uptake in c_uptakes:
        with model:
            if type=='ecGEM':
                model.reactions.get_by_id(modelParams['c_source']).bounds = 0, c_uptake
            else:
                model.reactions.get_by_id(modelParams['c_source']).bounds = -c_uptake, 0
            model.objective=modelParams['growth_id']
            growth=model.slim_optimize()
            growth_list.append(growth)
    growth_sim[type]=growth_list

df_result=pd.DataFrame(growth_sim,index=c_uptakes)

# load exp data
df_exp=pd.read_csv('examples/data/sce_growth_data.csv',sep='\t')


# plot the phenotype phase plane
import matplotlib.pyplot as plt
# set font style as Arial
plt.rcParams['font.sans-serif'] = 'Arial'
# set figure size
plt.figure(figsize=(5, 4))
# plot simulation result
plt.plot(df_result['gem'],label='Clasical GEM',linewidth=2)
plt.plot(growth_sim['ecGEM'],label='ecGEM',linewidth=2)
plt.plot(growth_sim['etfl'],label='EFL',linewidth=2)
# plot experiment result
plt.scatter(x=df_exp['qglucose'],y=df_exp['D(h−1)'],label='experimental data',marker='o',color='#434655')
plt.xlabel('Glucose uptake (mmol/gDW/h)',weight='bold',fontsize=16)
plt.ylabel('Specific growth rate (1/h)',weight='bold',fontsize=16)
# set 0 as
plt.xlim(0,25)
plt.ylim(0,1.4)
plt.legend()
# plt.title('Growth simulation of GEM,ecGEM and ETFL model')
plt.tight_layout()
plt.show()
