# -*- coding: utf-8 -*-
# date : 2024/4/17 
# author : wangh

'''How to evaluate the prediction results of strain design algorithms??
1. Calculate the experimental consistency of the prediction results by comparing the predicted gene targets with the experimental gene targets.
2.
'''

def calculate_exp_consistency(predict_result,exp_data,show=True):
    '''
    Calculate the experimental consistency of the prediction results by comparing the predicted gene targets with the experimental gene targets.
    '''
    predict_group=predict_result.groupby('actions')
    exp_group=exp_data.groupby('action')
    exp_consistency=dict()
    overall_exp_num=0
    overall_hit_num=0
    overall_predict_num=0
    for key in exp_group.groups.keys():
        exp_geneList=exp_group.get_group(key).index.tolist()
        predict_geneList=predict_group.get_group(key).index.tolist()
        hit_geneList=list(set(exp_geneList).intersection(set(predict_geneList)))
        overall_exp_num+=len(exp_geneList)
        overall_hit_num+=len(hit_geneList)
        overall_predict_num+=len(predict_geneList)
        exp_consistency[key]={'exp':exp_geneList,'predict':predict_geneList,'hit':hit_geneList,
                              'exp_num':len(exp_geneList),'hit_num':len(hit_geneList),
                              'consistency':len(set(exp_geneList).intersection(set(hit_geneList)))/len(exp_geneList)}
    exp_consistency['overall']={'exp_num':overall_exp_num,'hit_num':overall_hit_num,'predict_num':overall_predict_num,
                                'consistency':overall_hit_num/overall_exp_num}

    if show:
        # print(list(exp_consistency.keys()))
        for key in list(exp_consistency.keys()):
            print(f'{key}:')
            print(exp_consistency[key])


    return exp_consistency


import pandas as pd
from strainOptimizer.analysis.dataset import load_experiment_targets


# load result
result_ecmodel=pd.read_excel(r'examples/result/ecYeast_2-phenylethanol_gluc_10_ecFactory_result.xlsx',sheet_name='geneTable',index_col=0)
min_gene_result=result_ecmodel[result_ecmodel['minimal candidates set']==1]
result_group=min_gene_result.groupby('actions')
# load experiment data
exp_data=load_experiment_targets(product='2-phenylethanol')
# group according to action
exp_group=exp_data.groupby('action')

exp_consistency=dict()
overall_exp_num=0
overall_hit_num=0
overall_predict_num=0
for key in exp_group.groups.keys():
    exp_geneList=exp_group.get_group(key).index.tolist()
    predict_geneList=result_group.get_group(key).index.tolist()
    hit_geneList=list(set(exp_geneList).intersection(set(predict_geneList)))
    overall_exp_num+=len(exp_geneList)
    overall_hit_num+=len(hit_geneList)
    overall_predict_num+=len(predict_geneList)
    exp_consistency[key]={'exp':exp_geneList,'predict':predict_geneList,'hit':hit_geneList,
                          'exp_num':len(exp_geneList),'hit_num':len(hit_geneList),
                          'consistency':len(set(exp_geneList).intersection(set(hit_geneList)))/len(exp_geneList)}

exp_consistency['overall']={'exp_num':overall_exp_num,'hit_num':overall_hit_num,'predict_num':overall_predict_num,
                            'consistency':overall_hit_num/overall_exp_num}


