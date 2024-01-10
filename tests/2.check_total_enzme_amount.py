'''Check the total enzymes amount changing with different growth rates for ETFL model
'''
import pandas as pd
import numpy as np
from strainOptimizer.io import load_model
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes
from strainOptimizer.simulation.pprotFBA import pprotFBA_prot_conc

# load model
model=load_model('examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json',model_type='etfl',solver='optlang-gurobi')

# test_model=load_model('examples/models/ecoli/ecoli_core_curated.json',model_type='etfl',solver='optlang-gurobi')
# model=test_model

tol_enzs=0.1 # g/gDW
model=constrain_enzymes(model,tol_enzs,model_type='etfl')

c_source='r_1714'
c_uptakes=np.arange(1,11,1)
growth_id='r_2111'

all_enzIDlist=[enz.id for enz in model.enzymes]
exclusion = ['dummy_enzyme',  # 'rib', 'rib_mit', 'rnap', 'rnap_mit'
             ]
enzIDlist=[x for x in all_enzIDlist if x not in exclusion]

# fix glucose uptake rate
growths=[]
tol_enzs_amounts=[]
for c_uptake in c_uptakes:
    model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake
    gr=model.slim_optimize()
    growths.append(gr)
    all_enzs_concentration=pprotFBA_prot_conc(model=model,targetID=growth_id,c_source=c_source,c_uptake=c_uptake,enzymeIDlist=enzIDlist,model_type='etfl')
    total_enzymes=all_enzs_concentration.sum()
    tol_enzs_amounts.append(total_enzymes)
    print('max growth rate is %s and total proteins amount is %s when glucose uptake is %s'%(gr,total_enzymes,c_uptake))
df=pd.DataFrame({'growth':growths,'tol_enzs_amount':tol_enzs_amounts,'c_uptake':c_uptakes})

# plot growth and total enzymes amount with different glucose uptake rates
import matplotlib.pyplot as plt
fig,ax=plt.subplots()
ax.plot(df['c_uptake'],df['growth'],'r',label='max growth',marker='o')
ax.set_xlabel('glucose uptake rate')
ax.set_ylabel('growth rate(h-1)')
ax2=ax.twinx()
ax2.plot(df['c_uptake'],df['tol_enzs_amount'],'b',label='total enzymes amount',marker='o')
ax2.set_ylabel('total enzymes amount(g/gDW)')
# set legends for both ax and ax2 in the lower right corner
ax.legend(bbox_to_anchor=(1,0.2))
ax2.legend(bbox_to_anchor=(1,0.1))
# set title
plt.title('Fix total enzymes amount to %s g/gDW'%tol_enzs)
plt.tight_layout()
plt.show()


# compare with ecmodel
# load model
ecmodel=load_model('examples/models/yeast/ecYeastGEM_batch.xml',model_type='ecGEM')

model=ecmodel

# fix total enzymes amount
model=constrain_enzymes(model,tol_enzs,model_type='ecGEM')

c_source='r_1714_REV'
c_uptakes=np.arange(1,11,1)
growth_id='r_2111'

enzIDlist=[rxn.id for rxn in model.reactions if rxn.id.startswith('draw_prot_')]

growths=[]
tol_enzs_amounts=[]
for c_uptake in c_uptakes:
    with model:
        model.reactions.get_by_id(c_source).bounds=c_uptake,c_uptake
        gr=model.slim_optimize()
        growths.append(gr)
        # fix the growth rate
        model.reactions.get_by_id(growth_id).bounds=gr,gr
        prot_pool_rxnID='prot_pool_exchange'
        model.objective=prot_pool_rxnID
        model.objective_direction='min'
        total_enzymes=model.optimize().objective_value
        tol_enzs_amounts.append(total_enzymes)
    print('max growth rate is %s and total proteins amount is %s when glucose uptake is %s'%(gr,total_enzymes,c_uptake))
df_ec=pd.DataFrame({'growth':growths,'tol_enzs_amount':tol_enzs_amounts,'c_uptake':c_uptakes})

import matplotlib.pyplot as plt
fig,ax=plt.subplots()
ax.plot(df_ec['c_uptake'],df_ec['growth'],'r',label='max growth',marker='o')
ax.set_xlabel('glucose uptake rate')
ax.set_ylabel('growth rate(h-1)')
ax2=ax.twinx()
ax2.plot(df_ec['c_uptake'],df_ec['tol_enzs_amount'],'b',label='total enzymes amount',marker='o')
ax2.set_ylabel('total enzymes amount(g/gDW)')
# set legends for both ax and ax2 in the lower right corner
ax.legend(bbox_to_anchor=(1,0.2))
ax2.legend(bbox_to_anchor=(1,0.1))
# set title
plt.title('Fix total enzymes amount to %s g/gDW in ecGEM'%tol_enzs)
plt.tight_layout()
plt.show()