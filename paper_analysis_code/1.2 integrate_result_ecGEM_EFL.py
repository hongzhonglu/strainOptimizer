# -*- coding: utf-8 -*-
import pandas as pd
import os
import sys
sys.path.append(r'D:\code\github\strainOptimizer')
from strainOptimizer.analysis.dataset import calculate_exp_consistency,load_experiment_targets
from strainOptimizer.analysis.network import calculate_genetic_target_distance,count_adjacent_reactions_from_gene
from strainOptimizer.analysis.FCC import calculate_FCC_by_abundance
from strainOptimizer.io import load_model
import json
import pickle

# set work dir
os.chdir(r'D:\code\github\strainOptimizer')

# load
products_pathway_dict=json.load(open(r'data\products_pathway_targets.json'))

gem=load_model('examples/models/yeast/yeast-GEM.xml',model_type='gem')
with open('data/Sce_met_rxn_diGraph.pkl', "rb") as f:
    G = pickle.load(f)

# model=load_model('examples/models/yeast/yeast-GEM.xml',model_type='ecGEM')
productParam_dict={
    '2-phenylethanol':{'productName':'2-phenylethanol',
                       'targetID':'r_1589',
                       'ecGEM_filepath':'examples/models/yeast/ecYeastGEM_batch.xml'},
    'heme':{'productName':'heme',
                 'targetID':'EX_heme_a'
                 ,'ecGEM_filepath':'examples/models/yeast/heme_ecYeastGEM.xml'},
    'spermidine':{'productName':'spermidine',
                  'targetID':'r_2051',
                  'ecGEM_filepath':'examples/models/yeast/ecYeastGEM_batch.xml',
                  },
    'sclareol':{'productName':'sclareol',
                'targetID':'DM_sclareol_c',
                'ecGEM_filepath':'examples/models/yeast/sclareol_ecYeastGEM_batch.xml'},
    'ffa':{'productName':'free fatty acids',
                   'targetID':'r_2189',
                   'ecGEM_filepath':'examples/models/yeast/ecYeastGEM_batch.xml'},
}

full_result_dict = {'ecGEM': {}, 'etfl': {}}
slim_result_df = pd.DataFrame()
for file in os.listdir(r'analysis_code\results\ecGEM_vs_efl'):
    if file.endswith('.xlsx'):
        model_type = file.split('_')[0]
        product_name = file.split('_')[1]

        # # ignore gluc=1
        # if 'gluc_1_' in file:
        #     continue
        # load predicted result
        df = pd.read_excel(f'analysis_code/results/ecGEM_vs_efl/{file}', sheet_name='geneTable', index_col=0)

        # if too much targets, remove k score 0.3-0.5
        if df.shape[0] > 100:
            df = df[(df['k_score'] < 0.4) | (df['k_score'] > 1.05)]
            if df.shape[0] > 100:
                df = df[(df['k_score'] < 0.3) | (df['k_score'] > 1.1)]
                if df.shape[0] > 100:
                    df = df[(df['k_score'] < 0.2) | (df['k_score'] > 1.15)]
                    if df.shape[0] > 100:
                        df = df[(df['k_score'] < 0.1) | (df['k_score'] > 1.2)]

        # calculate FCCg, FCCp,adjacent reactions, distance to substrate and product
        model = load_model(productParam_dict[product_name]['ecGEM_filepath'], model_type='ecGEM')
        productID=productParam_dict[product_name]['targetID']
        for id in df.index:
            FCCg, FCCp = calculate_FCC_by_abundance(protID=id, model=model,
                                                    productID=productID,)
            substrate_distance, product_distance = calculate_genetic_target_distance(model=gem,
                                                                                     genetic_target=id,
                                                                                     substrate_ID='r_1714',
                                                                                     productID=productID)
            adjacent_reactions = count_adjacent_reactions_from_gene(model=gem, gene_id=id,G=G)
            df.loc[id, 'FCCg'] = FCCg
            df.loc[id, 'FCCp'] = FCCp
            df.loc[id, 'substrate_distance'] = substrate_distance
            df.loc[id, 'product_distance'] = product_distance
            df.loc[id, 'adjacent_reactions'] = adjacent_reactions


        # load experimental result
        if product_name == 'free fatty acids':
            product_name = 'ffa'
        exp_data = load_experiment_targets(product_name)

        # target_pathway_geneList=products_pathway_dict[product_name]
        # # remove target pathway from exp data and predict
        # exp_data=exp_data[~exp_data.index.isin(target_pathway_geneList)]
        # df=df[~df.index.isin(target_pathway_geneList)]

        level1_result = df
        level2_result = df[df['target_priority_leval'].isin([1, 2, 3])]
        level3_result = df[df['minimal candidates set'] == 1]

        l1_eval = calculate_exp_consistency(predict_result=level1_result, exp_data=exp_data, show=False)
        l2_eval = calculate_exp_consistency(predict_result=level2_result, exp_data=exp_data, show=False)
        l3_eval = calculate_exp_consistency(predict_result=level3_result, exp_data=exp_data, show=False)

        full_result_dict[model_type][product_name] = {'level1': l1_eval, 'level2': l2_eval, 'level3': l3_eval}

        # calculate inspiration index
        l1_consistency = l1_eval['overall']['consistency']
        l1_precision = l1_eval['overall']['precision']
        l1_predict_num = l1_eval['overall']['predict_num']
        l2_consistency = l2_eval['overall']['consistency']
        l2_precision = l2_eval['overall']['precision']
        l2_predict_num = l2_eval['overall']['predict_num']
        l3_consistency = l3_eval['overall']['consistency']
        l3_precision = l3_eval['overall']['precision']
        l3_predict_num = l3_eval['overall']['predict_num']

        # calculate exp consistency and accuracy
        # consistency_score=eval_result['overall']['consistency']
        # accuracy_score=eval_result['overall']['precision']
        # predict_num=eval_result['overall']['predict_num']
        # slim_result_dict[model_type][product_name]={'consistency':consistency_score,'accuracy':accuracy_score}
        # add new row to df
        new_row = pd.DataFrame(
            {'model_type': [model_type], 'product_name': [product_name], 'l1_consistency': [l1_consistency],
             'l1_precision': [l1_precision], 'l1_predict_num': [l1_predict_num], 'l2_consistency': [l2_consistency],
             'l2_precision': [l2_precision], 'l2_predict_num': [l2_predict_num], 'l3_consistency': [l3_consistency],
             'l3_precision': [l3_precision], 'l3_predict_num': [l3_predict_num]})
        slim_result_df = pd.concat([slim_result_df, new_row], ignore_index=True)