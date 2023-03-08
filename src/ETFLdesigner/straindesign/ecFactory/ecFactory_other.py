# -*- coding: utf-8 -*-
# date : 2023/3/7 
# author : wangh
import numpy as np
import pandas as pd
from etfl.io.json import load_json_model


def find_leaks(candidates, targetID, model):
    '''function to find reactions that consume the target.
    :param candidates: a pandas dataframe with the following columns:
        1. geneID: gene ID
        2. k_score: k-score of the gene
        3. actions: action type for the gene
    :param targetID: target product exchange reaction ID
    :param model: ETFL model
    :return: a pandas dataframe with the following columns:

    '''
    # assume target objective rxn must be the exchange rxn
    met_ex_rxn = model.reactions.get_by_id(targetID)
    met_e_id=met_ex_rxn.reactants[0].id
    met_c_id=met_e_id.replace('e','c')
    met_c=model.metabolites.get_by_id(met_c_id)
    with model:
        # set the target exchange reaction as objective
        model.objective = targetID
        # optimize
        sol = model.optimize()

    # Find reactions that consume the product in the cytoplasm
    met_c_sum = met_c.summary(solution=sol)
    consuming_rxnIDs= met_c_sum.consuming_flux.index.tolist()

    target_genes=[]
    for rxnID in consuming_rxnIDs:
        rxn=model.reactions.get_by_id(rxnID)
        # remove rxn to the product target
        met_idList=[met.id for met in list(rxn.metabolites.keys())]
        if met_e_id in met_idList:
            consuming_rxnIDs.remove(rxnID)
            continue
        else:
            # find gene taget for the rxn
            geneList=list(rxn.genes)
            geneIDlist=[gene.id for gene in geneList]
            # remove gene that has been included in candidates
            geneIDlist=[geneID for geneID in geneIDlist if geneID not in candidates.index.tolist()]

            target_genes=target_genes+geneIDlist

    df_new=pd.DataFrame({'geneID':target_genes,'k_score':np.nan,'actions':np.nan})
    df_new.set_index('geneID',inplace=True)
    df_new['k_score']=[0]*len(target_genes)
    df_new['actions']=['KO']*len(target_genes)
    print('%s leak rxn has been found and added to candidates.\n'%len(target_genes))
    # add to candidates
    candidates=candidates.append(df_new)

    return candidates


def remove_essential_targets(candidates,essential_path=r'code_etfl/ETFLdesigner/data/essential_genes.txt'):
    '''function to remove essential genes from candidates.
    :param candidates: a pandas dataframe with the following columns:
        1. geneID: gene ID
        2. k_score: k-score of the gene
        3. actions: action type for the gene
    :param essential_path: path to essential genes data(default:path for S.cerevisiae essential genes table)
    :return: a updated candidates dataframe
    '''
    # remove essential genes
    essentials = pd.read_csv(essential_path, sep='\t').Ids.str.strip()
    essential_targets = candidates.loc[candidates.index.isin(essentials)].index.tolist()
    print(f"Removing {len(essential_targets)} essential targets\n")
    candidates = candidates.drop(essential_targets)
    return candidates


def getMetGeneMatrix(model,geneIDlist=None):
    """
    Function that obtains a binary matrix in which rows represent metabolites and columns genes.
    Each non-zero coeffiecient represents a relationship between a gene and a metabolite
     (i.e. the metabolite is either consumed/produced by a reaction encoded by a given gene).
    Args:
    - model: (ETFL/ecGEM model) Model to obtain the matrix from
    - genes: (list) List of gene IDs or gene indexes to take from the model (Default: all genes in the model)

    Returns:
    - GeneMetMatrix: (numpy array) Boolean matrix representing the relationships metabolites and genes
    - metsConectivity: (pandas dataframe) Dataframe that indicates the amount of metabolites related to each gene
    - genesConectivity: (pandas dataframe) Dataframe that indicates the amount of genes related to each metabolite
    """
    # Manage exceptions
    if geneIDlist is None:
        genes = list(model.genes)
        geneIDlist = [g.id for g in genes]

    # build metGeneMatrix
    allmets_ID=[m.id for m in model.metabolites]
    metGeneMatrix=pd.DataFrame(index=allmets_ID,columns=geneIDlist)
    for gID in geneIDlist:
        g=model.genes.get_by_id(gID)
        rxnList=list(g.reactions)
        for rxn in rxnList:
            metList=list(rxn.metabolites.keys())
            for met in metList:
                metGeneMatrix.loc[met.id,gID]=1
    metGeneMatrix=metGeneMatrix.fillna(0)
    metGeneMatrix=metGeneMatrix.astype(int)
    # remove rows with all zeros
    metGeneMatrix=metGeneMatrix.loc[(metGeneMatrix!=0).any(axis=1)]

    # Calculate row and column sums in metGeneMatrix
    metsConectivity = pd.DataFrame(metGeneMatrix.sum(axis=1),columns=['gene_number'])
    genesConectivity = pd.DataFrame(metGeneMatrix.sum(axis=0),columns=['met_number'])

    return metGeneMatrix, metsConectivity, genesConectivity


def getGeneDepMatrix(metGeneMatrix,corr_threshold=0.99):
    """
    Function to estimate the linear dependencies between genes in a metabolites-genes binary matrix.
    Parameters:
    metGeneMatrix (pd.Dataframe): Binary matrix indicating the relationships,rows represent metabolites and columns represent genes.
    corr_threshold (float): Correlation threshold to consider two genes as linearly dependent (Default: 0.99).

    Returns:
    gene_equal_matrix (pd.Dataframe): A matrix indicating the linear dependencies between genes.1 means the two genes are equal that they connetct to the same metabolites.
    independat_genes (pd.Series): a series indicate the independent genes (1),which are not equal to any other genes.

    """

    # get all target gene IDs
    geneList = metGeneMatrix.columns.tolist()

    # build a correlation matrix: calculate pearson correlation between genes based on related metabolites
    corr_matrix = np.corrcoef(metGeneMatrix.T)
    df_gene_corr_matrix = pd.DataFrame(corr_matrix, index=geneList, columns=geneList)
    gene_equal_matrix=df_gene_corr_matrix[df_gene_corr_matrix>corr_threshold]
    # fill nan as 0
    gene_equal_matrix=gene_equal_matrix.fillna(0)

    # estimate weather the genes are linearly dependent or not based on the correlation threshold.
    dep_genes = np.abs(corr_matrix) > corr_threshold

    #get the independent genes
    # remove the correlationship with the gene itself : Set the diagonal to False
    np.fill_diagonal(dep_genes, False)
    # identify the independent genes which has no linear relationship with any other genes
    indep_genes = np.sum(dep_genes, axis=0) == 0
    independant_genes = pd.Series(indep_genes, index=geneList)
    # fill True as 1 and False as 0
    independant_genes = independant_genes.replace({True: 1, False: 0})

    # return independant_genes, gene_equal_matrix, df_gene_corr_matrix #check the result by return the correlation matrix simultaneously.
    return independant_genes, gene_equal_matrix


def getGenesGroups(gene_equal_Matrix):
    """
    Function to get the groups of genes that are linearly dependent.
    Parameters:
    gene_equal_matrix (pd.Dataframe):
        A matrix indicating the linear dependencies between genes.1 means the two genes are equal that they connetct to the same metabolites.
    returns:
    geneGroups (pd.Series): target genes groups which have the same related metabolites.
    """
    # get the gene groups which have the same equal genes
    geneGroups0 = gene_equal_Matrix.groupby(gene_equal_Matrix.columns.tolist(), axis=1).groups
    # convert the gene groups to list
    geneGroups0 = {k: list(v) for k, v in geneGroups0.items()}
    i=1
    geneGroups= {}
    gene_numbers=[]
    for k, v in geneGroups0.items():
        genes_numb=len(v)
        if genes_numb>1:
            geneGroups[i]=v
            gene_numbers.append(genes_numb)
            i+=1

    geneGroups=pd.Series(geneGroups)

    return geneGroups


if __name__ == '__main__':
    test_model= load_json_model('models/ecoli_core.json')
    # yefl=load_json_model('models/yeast8_cEFL_2584_enz_128_bins__20221115_120238.json')
    model=test_model
    geneTable=pd.DataFrame(columns=['geneID','k_score','actions'])
    # geneTable=pd.read_excel('code_etfl/ETFLdesigner/output/yefl_2PE_design_results.xlsx',sheet_name='geneTable')
    targetID='r_1589'
    candidates1=find_leaks(candidates=geneTable, targetID=targetID, model=model)

    # test getMetGeneMatrix
    metGeneMatrix, metsConectivity, genesConectivity = getMetGeneMatrix(model)
    # remove all 0 rows and columns
    metGeneMatrix = metGeneMatrix.loc[(metGeneMatrix != 0).any(axis=1)]
    metGeneMatrix = metGeneMatrix.loc[:, (metGeneMatrix != 0).any(axis=0)]

    # test getGeneDepMatrix
    for i in np.linspace(0.1,0.95,20):
        independant_genes, gene_equal_matrix= getGeneDepMatrix(metGeneMatrix,corr_threshold=i)
        print(f'corr_threshold={i},independant_genes={independant_genes.sum()}')

    # test getGenesGroups
    geneGroups=getGenesGroups(gene_equal_matrix)