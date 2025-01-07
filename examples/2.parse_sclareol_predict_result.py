# -*- coding: utf-8 -*-
# date : 2024/12/24 
import pandas as pd
from strainOptimizer.analysis import dataset

def parse_hit_geneName(consistency_result,experimental_data):
    kd_hit_geneList=consistency_result['KD']['hit']
    oe_hit_geneList=consistency_result['OE']['hit']
    print(f'Up regulation experimental hit genes: {len(oe_hit_geneList)}')
    print(f"Down regulation experimental hit gene: {len(kd_hit_geneList)}")
    hit_geneList=kd_hit_geneList+oe_hit_geneList
    hit_data=experimental_data.loc[hit_geneList]
    return hit_data

# load ecFactory result
df_ecYeast=pd.read_excel(r'examples/result/ecGEM_sclareol_gluc_1_ecFactory_result.xlsx',sheet_name='geneTable',index_col=0)
df_efl=pd.read_excel('examples/result/etfl_sclareol_gluc_5_ecFactory_result.xlsx',sheet_name='geneTable',index_col=0)

df_ecYeast['geneName']=dataset.gene_id_to_name(df_ecYeast.index)
df_efl['geneName']=dataset.gene_id_to_name(df_efl.index)

# load experimental dataset
exp_data=dataset.load_experiment_targets(product='sclareol')

# Calculate experimental consistency
ecYeast_consistency=dataset.calculate_exp_consistency(predict_result=df_ecYeast,exp_data=exp_data)
efl_consistency=dataset.calculate_exp_consistency(predict_result=df_efl,exp_data=exp_data)

ecYeast_hit_gene=parse_hit_geneName(consistency_result=ecYeast_consistency,experimental_data=exp_data)
efl_hit_gene=parse_hit_geneName(consistency_result=efl_consistency,experimental_data=exp_data)


