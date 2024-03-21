# -*- coding: utf-8 -*-
# date : 2024/3/20 
# author : wangh
'''Test MOPA simulation for ETFL/ecGEM model.
'''
import pandas as pd
from strainOptimizer.simulation import mopa,moma
from strainOptimizer.io import load_model
from strainOptimizer.simulation.pprotFBA import ppFBA2

solver='optlang-gurobi'

# test ecGEM
ecYeast = load_model("examples/models/yeast/heme_ecYeastGEM.xml",model_type='ecGEM',solver=solver)
model=ecYeast

model.tolerance = 1e-9
product_name = 'heme_a'
product_id = 'EX_heme_a'
growth = 'r_2111'
c_source = 'r_1714_REV'
c_uptake = 10
prot_pool='prot_pool_exchange'


# set wild type condition as the reference
with model:
    ref_sol=ppFBA2(model=model, targetID=growth, c_source=c_source, c_uptake=c_uptake, model_type='ecGEM')


with model:
    model.objective = product_id
    product=model.slim_optimize()*0.5
mutant = model.copy()
mutant.reactions.get_by_id(product_id).bounds = (product, float('inf'))

# test MOPA
sol = mopa(model=mutant, reference_solution=ref_sol, linear=False,model_type='ecGEM')
sol_linear=mopa(model=mutant, reference_solution=ref_sol, linear=True,model_type='ecGEM')

sol_moma =  moma(model=mutant, reference_solution=ref_sol, linear=False,model_type='ecGEM')
sol_moma_linear = moma(model=mutant, reference_solution=ref_sol, linear=True,model_type='ecGEM')

df_compare=pd.DataFrame({'reference':ref_sol.fluxes,'MOPA':sol.fluxes,'MOPA_linear':sol_linear.fluxes,'MOMA':sol_moma.fluxes,'MOMA_linear':sol_moma_linear.fluxes})
# plot the scatter plot
import matplotlib.pyplot as plt
fig,axes=plt.subplots(2,3,figsize=(15,10))
count=0
for i in df_compare.columns:
    # plt.scatter(df_compare['reference'],df_compare[i],label=i,ax=ax[count//3,count%3])
    axes[count//3,count%3].scatter(df_compare['reference'],df_compare[i],label=i)
    axes[count//3,count%3].set_title(i,fontweight='bold')
    axes[count//3,count%3].set_xlabel('reference')
    axes[count//3,count%3].set_ylabel(i)
    count+=1
# set title
fig.suptitle('MOPA and MOMA simulation for ecGEM',fontweight='bold',fontsize=20)
fig.show()


# test ETFL modle
etfl=load_model("examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json", solver=solver,model_type='etfl')
model=etfl

# set tolerance as 1e-9
model.tolerance = 1e-9

# productID='EX_succ_e'
# growth_id='Biomass_Ecoli_core'
# glucose_id='EX_glc__D_e'
# dummy_enzID='dummy_enzyme'

# for ETFL model
productID = 'r_1589'   # 2-PE production rxn
growth_id = 'r_2111'
c_source = 'r_1714'
c_uptake = 10
dummy_enzID = 'dummy_enzyme'


# 1. get reference solution
ref_sol=ppFBA2(model=model, targetID=growth_id, c_source=c_source, c_uptake=c_uptake, model_type='etfl')


with model:
    model.objective = productID
    product=model.slim_optimize()*0.5
mutant = model.copy()
mutant.reactions.get_by_id(productID).bounds = (product, 1000)
mutant.reactions.get_by_id(c_source).bounds = -c_uptake, 0
mutant.reactions.get_by_id(growth_id).bounds = 0, 1000
mutant.slim_optimize()
# test MOPA
sol = mopa(model=mutant, reference_solution=ref_sol, linear=False,model_type='etfl')
sol_linear=mopa(model=mutant, reference_solution=ref_sol, linear=True,model_type='etfl')

sol_moma =  moma(model=mutant, reference_solution=ref_sol, linear=False,model_type='etfl')
sol_moma_linear = moma(model=mutant, reference_solution=ref_sol, linear=True,model_type='etfl')

df_compare=pd.DataFrame({'reference':ref_sol.fluxes,'MOPA':sol.fluxes,'MOPA_linear':sol_linear.fluxes,'MOMA':sol_moma.fluxes,'MOMA_linear':sol_moma_linear.fluxes})
df_compare['modified_MOPA_linear']=df_compare['MOPA_linear'].apply(lambda x: x if abs(x)<40 else 40)
# plot the scatter plot
import matplotlib.pyplot as plt
fig,axes=plt.subplots(2,3,figsize=(15,10))
count=0
for i in df_compare.columns:
    # plt.scatter(df_compare['reference'],df_compare[i],label=i,ax=ax[count//3,count%3])
    axes[count//3,count%3].scatter(df_compare['reference'],df_compare[i],label=i)
    axes[count//3,count%3].set_title(i,fontweight='bold')
    axes[count//3,count%3].set_xlabel('reference')
    axes[count//3,count%3].set_ylabel(i)
    count+=1
# set title
fig.suptitle('MOPA and MOMA simulation for ETFL',fontweight='bold',fontsize=20)
fig.show()
