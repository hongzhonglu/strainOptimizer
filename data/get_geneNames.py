# -*- coding: utf-8 -*-
import pandas as pd

df_annotation=pd.read_excel(r'data/s288c_gene_annotation_score.xlsx',index_col=0)

df_geneName_id=df_annotation[['Gene names','geneID']]
# remove geneID NAN
df_geneName_id=df_geneName_id.dropna(axis=0,how='any')
df_geneName_id['geneNames']=df_geneName_id['Gene names'].apply(lambda x: ';'.join(str(x).split(' ')))
df_geneName_id['geneName']=df_geneName_id['Gene names'].apply(lambda x: str(x).split(' ')[0])

# drop columns
df_geneName_id=df_geneName_id.drop(columns=['Gene names'])

df_geneName_id.reset_index(inplace=True)
# set geneID as index
df_geneName_id.set_index('geneID',inplace=True)

# save result
df_geneName_id.to_csv(r'data/s288c_geneNames.csv')