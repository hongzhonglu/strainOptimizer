# -*- coding: utf-8 -*-
import pandas as pd
import os
import sys
sys.path.append(r'D:\code\github\strainOptimizer')
from strainOptimizer.analysis.dataset import calculate_exp_consistency,load_experiment_targets
from strainOptimizer.analysis.network import calculate_genetic_target_distance,count_adjacent_reactions_from_gene
from strainOptimizer.analysis.FCC import calculate_FCC_by_abundance,calculate_FCC_by_kcat
from strainOptimizer.io import load_model
import json
import pickle
import tqdm

def supply_analysis(df_geneTable,productParam_dict,gem):
    '''Supply analysis for the given gene table and ecGEM model.
    Args:
        df_geneTable (pd.DataFrame): The gene table with target genes.
        ecGEM_filepath (str): The file path of the ecGEM model.
        productName (str): The name of the product.
        targetID (str): The target reaction ID for the product.
        model_type (str): The type of the model, either 'ecGEM' or 'etfl'.
    Returns:
        pd.DataFrame: A DataFrame containing the supply analysis results.
    '''

    # control the total target number
    # remove k_score < 0
    df_geneTable = df_geneTable[df_geneTable['k_score'] >= 0].copy()

    if df_geneTable.shape[0] > 100:
        df_geneTable = df_geneTable[(df_geneTable['k_score'] < 0.4) | (df_geneTable['k_score'] > 1.05)].copy()
        if df_geneTable.shape[0] > 100:
            df_geneTable = df_geneTable[(df_geneTable['k_score'] < 0.3) | (df_geneTable['k_score'] > 1.1)].copy()
            if df_geneTable.shape[0] > 100:
                df_geneTable = df_geneTable[(df_geneTable['k_score'] < 0.2) | (df_geneTable['k_score'] > 1.15)].copy()
                if df_geneTable.shape[0] > 100:
                    df_geneTable = df_geneTable[(df_geneTable['k_score'] < 0.1) | (df_geneTable['k_score'] > 1.2)].copy()
                    if df_geneTable.shape[0] > 100:
                        # only keep top50 and last50 target, and remove duplicate
                        df_geneTable = df_geneTable.sort_values(by='k_score',ascending=False)
                        top50_list=df_geneTable.index[:45].tolist()
                        last50_list=df_geneTable.index[-45:].tolist()
                        selected_list=list(set(top50_list+last50_list))
                        df_geneTable = df_geneTable[df_geneTable.index.isin(selected_list)]
    
    # calculate FCCg, FCCp
    model=load_model(productParam_dict['ecGEM_filepath'],model_type='ecGEM')
    for id in df_geneTable.index:
        try:
            gene = model.genes.get_by_id(id)
        except:
            continue
        protID = None
        for rxn in gene.reactions:
            if 'draw_prot' in rxn.id:
                protID = rxn.id
                break
        if protID is not None:
            FCCg,FCCp=calculate_FCC_by_abundance(protID=protID,model=model,productID=productParam_dict['targetID'],delta_conc=1)
            # FCCg,FCCp=calculate_FCC_by_kcat(protID=protID,model=model,productID=productParam_dict['targetID'],delta_kcat=1)
            print(f'Gene {id}: FCCg={FCCg}, FCCp={FCCp}')
            df_geneTable.loc[id,'FCCg']=FCCg
            df_geneTable.loc[id,'FCCp']=FCCp

        try:
            go=gem.genes.get_by_id(id)
        except:
            print(f'Gene {id} not found in the GEM model.')
            continue
        # calculate distance to substrate and product
        distance_to_substrate,distance_to_product=calculate_genetic_target_distance(model=gem,genetic_target=id,substrate_ID='r_1714',productID=productParam_dict['targetID'])
        df_geneTable.loc[id,'distance_to_substrate']=distance_to_substrate
        df_geneTable.loc[id,'distance_to_product']=distance_to_product
        df_geneTable.loc[id,'distance']=min(distance_to_substrate,distance_to_product)

        # calculate adjacent reactions
        count,adj_rxnList=count_adjacent_reactions_from_gene(model=gem,gene_id=id)
        df_geneTable.loc[id,'count']=count
        df_geneTable.loc[id,'adj_rxnList']=str(adj_rxnList)
    return df_geneTable

# set work dir
os.chdir(r'D:\code\github\strainOptimizer')

gem=load_model('examples/models/yeast/yeast-GEM.xml',model_type='gem')
with open('data/slim_Sce_met_rxn_diGraph.pkl', "rb") as f:
    G = pickle.load(f)

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

file_dir=r'analysis_code/results/moma_mopa_pfba'
result_dict={}
for file in tqdm.tqdm(os.listdir(file_dir)):
    if file.endswith('.xlsx'):
        product_name=file.split('_')[1]

    # load experimental result
    if product_name=='free fatty acids':
        product_name='ffa'

    productParam=productParam_dict[product_name]

    df = pd.read_excel(f'{file_dir}/{file}', sheet_name='geneTable', index_col=0)

    if 'FCCg' in df.columns:
        continue

    df_add= supply_analysis(df_geneTable=df, productParam_dict=productParam,gem=gem)
    result_dict[product_name]=df_add

    # save to excel
    with pd.ExcelWriter(f'{file_dir}/{file}', engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
        df_add.to_excel(writer, sheet_name='geneTable')



