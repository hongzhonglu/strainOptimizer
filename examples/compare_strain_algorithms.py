# -*- coding: utf-8 -*-
# date : 2024/4/17 
# author : wangh
# file : compare_algorithm.py
# project : strainOptimizer
'''Compare different algorithms for strain design.'''
from strainOptimizer.analysis.dataset import calculate_exp_consistency,load_experiment_targets
import pandas as pd

def run_evaluation_process(predict_results_dict,exp_data):
    for key in predict_results_dict.keys():
        print(f'{key}:')
        result=calculate_exp_consistency(predict_results_dict[key],exp_data,show=False)
        print(result['overall'])
        print('---------------------')


# 2-phenylethanol
product_name='2-phenylethanol'
# load experiment data
exp_data_2pe=load_experiment_targets(product_name)

# load ecFactory result
ecfactory_candidates_l3=pd.read_csv('reference/ecFactory-main/tutorials/results/2_phenylethanol_targets/candidates_L3.txt',sep='\t',index_col=0)
ecfactory_candidates_l1=pd.read_csv('reference/ecFactory-main/tutorials/results/2_phenylethanol_targets/candidates_L1.txt',sep='\t',index_col=0)
ecfactory_candidates_dict={'l1':ecfactory_candidates_l1,
    'l3':ecfactory_candidates_l3}


ecmodel_result=pd.read_excel(r'examples/result/ecYeast_2-phenylethanol_gluc_1_ecFactory_result.xlsx',sheet_name='geneTable',index_col=0)
ecmodel_candidates_l3=ecmodel_result[ecmodel_result['minimal candidates set']==1]
ecmodel_candidates_l2=ecmodel_result[ecmodel_result['target_priority_leval'].isin([1,2,3])]
ecmodel_candidates_l1=ecmodel_result
ecmodel_result_dict={'l1':ecmodel_candidates_l1,
    'l2':ecmodel_candidates_l2,
                     'l3':ecmodel_candidates_l3}

etfl_result=pd.read_excel(r'examples/result/yefl_2-phenylethanol_gluc_1_ecFactory_result.xlsx',sheet_name='geneTable',index_col=0)
etfl_candidates_l3=etfl_result[etfl_result['minimal candidates set']==1]
etfl_candidates_l2=etfl_result[etfl_result['target_priority_leval'].isin([1,2,3])]
etfl_candidates_l1=etfl_result
etfl_result_dict={'l1':etfl_candidates_l1,
    'l2':etfl_candidates_l2,
                     'l3':etfl_candidates_l3}


print('ecFactory:')
run_evaluation_process(ecfactory_candidates_dict,exp_data_2pe)
print('etfl:')
run_evaluation_process(etfl_result_dict,exp_data_2pe)
print('ecmodel:')
run_evaluation_process(ecmodel_result_dict,exp_data_2pe)

# for heme production
product_name='heme'
# load experiment data
exp_data_heme=load_experiment_targets(product_name)

# load ecFactory result
ecfactory_candidates_l3=pd.read_csv(r'reference/ecFactory-main/tutorials/results/heme_targets/candidates_L3.txt',sep='\t',index_col=0)
ecfactory_candidates_l1=pd.read_csv(r'reference/ecFactory-main/tutorials/results/heme_targets/candidates_L1.txt',sep='\t',index_col=0)
ecfactory_candidates_dict={'l1':ecfactory_candidates_l1,
    'l3':ecfactory_candidates_l3}

# load ecmodel result
ecmodel_result=pd.read_excel(r'examples/result/ecYeast_heme_a_gluc_1_ecFactory_result.xlsx',sheet_name='geneTable',index_col=0)
ecmodel_candidates_l3=ecmodel_result[ecmodel_result['minimal candidates set']==1]
ecmodel_candidates_l2=ecmodel_result[ecmodel_result['target_priority_leval'].isin([1,2,3])]
ecmodel_candidates_l1=ecmodel_result
ecmodel_result_dict={'l1':ecmodel_candidates_l1,
    'l2':ecmodel_candidates_l2,
                     'l3':ecmodel_candidates_l3}

print('ecmodel:')
run_evaluation_process(ecmodel_result_dict,exp_data_heme)
print('ecFactory:')
run_evaluation_process(ecfactory_candidates_dict,exp_data_heme)



