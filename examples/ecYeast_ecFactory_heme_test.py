# -*- coding: utf-8 -*-
# date : 2023/4/6 
# author : wangh
# file : ecYeast_ecFactory_heme_test.py
# project : etfl
import pandas as pd
import numpy as np
from ETFLdesigner.io.ecModel import load_ecmodel
from ETFLdesigner.strain_design.ecFactory import run_ecFactory

# load heme a product model
heme_ecYeast=load_ecmodel('models/heme_ecYeastGEM.xml')
product_name='heme a'
product_id='EX_heme_a'
model=heme_ecYeast

# calculate the max biomass yield
# set tolerance
model.tolerance=1e-9

# parameters for ecYeast
c_source="r_1714_REV"      # glucose exchange rxn
growth_id='r_2111'
model.objective=growth_id# biomass rxn
c_uptake=1
model.reactions.get_by_id('r_1714_REV').bounds=1,1
gluc_MW=0.180156  # g/mmol
max_yield=model.slim_optimize()/(c_uptake*gluc_MW) # gDW / gGluc
expYield=0.122
alphaLims=(0.5*expYield,2*expYield)

# prepare parameters for ecFactory algorithm
modelParam=pd.Series()
modelParam['targetID']=product_id
modelParam['c_source']=c_source
modelParam['c_uptake']=c_uptake
modelParam['productName']=product_name
action_thresholds=[0.05,0.3,1.05]     # rules for overexpression, knockout and knockdown

# set a timer
import time
start_time = time.time()
print('start time:',start_time)
# run ecFSEOF
results=run_ecFactory.run_ecFactory_design(model=model, modelParam=modelParam, expYield=expYield,alphaLims=alphaLims,action_thresholds=action_thresholds,remove_essential=True,model_type='ecGEM')
end_time = time.time()
print('end time:',end_time)
print('time cost:',end_time-start_time)

# print FSEOF results
print('FSEOF results:')
genetable=results['geneTable']
print('OE:',genetable[genetable['actions']=='OE'].shape[0])
print('KD:',genetable[genetable['actions']=='KD'].shape[0])
print('KO:',genetable[genetable['actions']=='KO'].shape[0])

# EUVA results
print('EUVA results:')
print('leval 1 candidates:',genetable[genetable['target_priority_leval']==1].shape[0])
print('leval 2 candidates:',genetable[genetable['target_priority_leval']==2].shape[0])
print('leval 3 candidates:',genetable[genetable['target_priority_leval']==3].shape[0])


# minimal combined set
min_result=results['min_set_analysis_result']
min_set_list=min_result[min_result['score']<0.999].index.tolist()
fseof_list=results['geneTable'].index.tolist()
euva_list=results['geneTable'][results['geneTable']['target_priority_leval'].isin([1,2,3])].index.tolist()


# compare result with ecFactory
# load ecFactory predicted result
ecfactory_candidates_l3=pd.read_csv('code_etfl/reference/ecFactory/results/heme_targets/candidates_L3.txt',sep='\t',index_col=0).index.tolist()
ecfactory_candidates_l1=pd.read_csv('code_etfl/reference/ecFactory/results/heme_targets/candidates_L1.txt',sep='\t',index_col=0).index.tolist()

common_l3=set(min_set_list).intersection(set(ecfactory_candidates_l3))
comon_l1=set(min_set_list).intersection(set(ecfactory_candidates_l1))
fseof_common_l3=set(fseof_list).intersection(set(ecfactory_candidates_l3))
fseof_common_l1=set(fseof_list).intersection(set(ecfactory_candidates_l1))
euva_common_l3=set(euva_list).intersection(set(ecfactory_candidates_l3))
euva_common_l1=set(euva_list).intersection(set(ecfactory_candidates_l1))
print('--------------------------compare result with ecFactory result')
print('minimal sets vs l3:',len(common_l3))
print('minimal sets vs l1:',len(comon_l1))
print('fseof vs l3:',len(fseof_common_l3))
print('fseof vs l1:',len(fseof_common_l1))
print('euva vs l3:',len(euva_common_l3))
print('euva vs l1:',len(euva_common_l1))


# compare with experimental result
print('--------------------------compare predicted candidates with experimental data')
# load heme experimental data
heme_exp=pd.read_excel('ETFLdesigner/data/heme_experimental_data.xlsx',index_col=0)
exp_list=heme_exp.index.tolist()
common_exp=set(min_set_list).intersection(set(exp_list))
all_common_exp=set(fseof_list).intersection(set(exp_list))
euva_common_exp=set(euva_list).intersection(set(exp_list))
print('minimal sets vs exp:',len(common_exp))
print('fseof vs exp:',len(all_common_exp))
print('euva vs exp:',len(euva_common_exp))



# save results into excel file
for key in results.keys():
    if isinstance(results[key],np.ndarray):
        results[key]=pd.Series(results[key])
    # 检查是否为dict
    if isinstance(results[key],dict):
        results[key]=pd.Series(results[key])
# with pd.ExcelWriter(f'code_etfl/ETFLdesigner/output/ecYeast_{product_name}_gluc_{c_uptake}_ecFactory_result.xlsx') as writer:
with pd.ExcelWriter(f'code_etfl/ETFLdesigner/output/ecYeast_{product_name}_gluc_{c_uptake}_ecFSEOF_result.xlsx') as writer:
    for key in results.keys():
        results[key].to_excel(writer, sheet_name=key)


# save geneTable
genetable=results['geneTable']
genetable.to_excel(f'ETFLdesigner/examples/siwei_heme_FSEOF/why_heme_geneTable.xlsx')
results['v_matrix'].to_excel(f'ETFLdesigner/examples/siwei_heme_FSEOF/why_heme_v_matrix.xlsx')
results['k_matrix'].to_excel(f'ETFLdesigner/examples/siwei_heme_FSEOF/why_heme_k_matrix.xlsx')


