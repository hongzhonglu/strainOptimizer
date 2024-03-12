# -*- coding: utf-8 -*-
import pandas as pd

def genelist_to_enzymelist(model,genelist,model_type='etfl'):
    '''
    get the enzyme list from a gene list
    para:
        model: must be ETFL model/ ecGEM model
        genelist: a list of gene ID
    return:
        a list of enzyme ID
    '''
    if model_type=='etfl':
        # get all enzyme to gene dict
        all_enzIDlist=[enz.id for enz in model.enzymes]
        enz_geneDict={}
        for enzID in all_enzIDlist:
            enz_i_genelist=list(model.enzymes.get_by_id(enzID).composition.keys())
            enz_i_gene_list_to_str='|'.join(enz_i_genelist)
            enz_geneDict[enzID]=enz_i_gene_list_to_str
        df_enz_gene=pd.Series(enz_geneDict)
        # find target enzymes list
        enzlist=[]
        gene_enz_dict={}
        for gene in genelist:
            enzymes=df_enz_gene[df_enz_gene.str.contains(gene)].index.tolist()
            gene_enz_dict[gene]=enzymes
            enzlist=enzlist+enzymes
        enzlist=list(set(enzlist))

    elif model_type=='ecGEM':
        # gene_to_prot_dict
        gene_to_prot_dict = {}
        for g in model.genes:
            geneID = g.id
            for rxn in g.reactions:
                if 'draw_prot_' in rxn.id:
                    draw_prot_rxnID = rxn.id
                    gene_to_prot_dict[geneID] = draw_prot_rxnID
        all_enz_genes= list(gene_to_prot_dict.keys())
        # find target enzymes list
        enzlist=[]
        gene_enz_dict={}
        for geneID in genelist:
            if geneID not in all_enz_genes:
                gene_enz_dict[geneID]='no enzyme'
            else:
                enzID=gene_to_prot_dict[geneID]
                enzlist.append(enzID)
                gene_enz_dict[geneID]=[enzID]

    return enzlist,gene_enz_dict