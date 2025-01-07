# -*- coding: utf-8 -*-
'''Load standard datasets foe strain design algorithm evaluation
'''

import pandas as pd
import os

def load_experiment_targets(product:str, data_dir='data/experiment_targets'):
    '''Load experiment targets for a specific product
    '''
    available_products = [f.replace('_exp_targets.tsv','') for f in os.listdir(data_dir) if f.endswith('_exp_targets.tsv')]
    if product not in available_products:
        print('Available products:', available_products)
        raise ValueError('The product %s is not available!' % product)
    else:
        df = pd.read_csv(os.path.join(data_dir, product+'_exp_targets.tsv'), sep='\t', index_col=0)
        return df


def calculate_exp_consistency(predict_result, exp_data, show=True):
    '''
    Calculate the experimental consistency of the prediction results by comparing the predicted gene targets with the experimental gene targets.
    '''
    predict_result= predict_result[predict_result['action'].isin(['OE', 'KD', 'KO'])]
    predict_group = predict_result.groupby('action')
    exp_group = exp_data.groupby('action')
    exp_consistency = dict()
    overall_exp_num = 0
    overall_hit_num = 0
    overall_predict_num = 0
    for key in exp_group.groups.keys():
        exp_geneList = exp_group.get_group(key).index.tolist()
        try:
            predict_geneList = predict_group.get_group(key).index.tolist()
        except:
            predict_geneList = []
        hit_geneList = list(set(exp_geneList).intersection(set(predict_geneList)))
        overall_exp_num += len(exp_geneList)
        overall_hit_num += len(hit_geneList)
        overall_predict_num += len(predict_geneList)
        exp_consistency[key] = {'exp': exp_geneList, 'predict': predict_geneList, 'hit': hit_geneList,
                                'exp_num': len(exp_geneList), 'hit_num': len(hit_geneList),
                                'consistency': len(set(exp_geneList).intersection(set(hit_geneList))) / len(
                                    exp_geneList)}
    exp_consistency['overall'] = {'exp_num': overall_exp_num, 'hit_num': overall_hit_num,
                                  'predict_num': overall_predict_num,
                                  'consistency': overall_hit_num / overall_exp_num}

    if show==True:
        for key in exp_consistency.keys():
            print(f'{key}:')
            print(exp_consistency[key])

    return exp_consistency


def gene_id_to_name(geneIDlist,annotation_file=r'data/s288c_geneNames.csv'):
    df=pd.read_csv(annotation_file,index_col=0)
    df_geneName=df[df.index.isin(geneIDlist)]['geneName']
    return df_geneName