# -*- coding: utf-8 -*-
# date : 2024/5/22 
# author : wangh
# file : build_lexicon_for_ecmodel.py
# project : strainOptimizer
from cobra.io import read_sbml_model,load_matlab_model
import pandas as pd

curated_model=load_matlab_model(r'D:\code\github\etfl\yetfl\info_yeast\Y8_3_4_mod_curated.mat')
ecModel=read_sbml_model('examples/models/yeast/ecYeastGEM_batch.xml')

# load curated_model lexicon
df_lexcicon_curated = pd.read_csv(r'data/yeast_thermo/yeast_lexicon.csv',index_col=0)

gem_met1=curated_model.metabolites.get_by_id('s_0001_c')
ec_met1=ecModel.metabolites.get_by_id('s_0001[ce]')

df_lexcicon_curated['name']=df_lexcicon_curated.index.map(lambda x: curated_model.metabolites.get_by_id(x).name)
# remove the [xxx] in the end of string in name column
df_lexcicon_curated['name']=df_lexcicon_curated['name'].map(lambda x: x.split('[')[0].strip())
# group seed_id according to name
seed_met=df_lexcicon_curated.groupby('name')['seed_id'].apply(lambda x: x.tolist()[0])

# build ecGEM lexicon
lexicon_dict=dict()
for met in ecModel.metabolites:
    id=met.id
    name=met.name.split('[')[0].strip()
    if name in seed_met.index:
        lexicon_dict[id]=seed_met[name]
df_lexicon_ecmodel=pd.DataFrame.from_dict(lexicon_dict,orient='index',columns=['seed_id'])

# save result
df_lexicon_ecmodel.to_csv(r'data/yeast_thermo/ecYeast_lexicon.csv')


# build lexicon for yeast9
yeast9=read_sbml_model('examples/models/yeast/yeast-GEM.xml')

df_lexicon_ecmodel=pd.read_csv(r'data/yeast_thermo/ecYeast_lexicon.csv',index_col=0)
df_lexcicon_yeast9=df_lexicon_ecmodel.copy()
df_lexcicon_yeast9['id']=df_lexcicon_yeast9.index.map(lambda x: x.split('[')[0].strip())
# set id as index
df_lexcicon_yeast9.set_index('id',inplace=True)
metIDlist=[x.id for x in yeast9.metabolites]

# remove metabolites not in yeast9
df_lexcicon_yeast9=df_lexcicon_yeast9.loc[df_lexcicon_yeast9.index.isin(metIDlist)]

# add new metabolites id
to_add_dict={}
for id in metIDlist:
    if id not in df_lexcicon_yeast9.index:
        met=yeast9.metabolites.get_by_id(id)
        if 'kegg.compound' in met.annotation.keys():
            to_add_dict[id]=met.annotation['kegg.compound']
        elif 'chebi' in met.annotation.keys():
            to_add_dict[id]=met.annotation['chebi']
        else:
            to_add_dict[id]=met.name
df_to_add=pd.DataFrame.from_dict(to_add_dict,orient='index',columns=['seed_id'])
df_lexcicon_yeast9=pd.concat([df_lexcicon_yeast9,df_to_add])

# save result
df_lexcicon_yeast9.to_csv(r'data/yeast_thermo/yeast9_lexicon.csv')