# -*- coding: utf-8 -*-
# date : 2023/2/6 
# author : wangh
# file : concentration_predict.py
# project : etfl
# predict enzyme concentration ,mRNA concentration
from etfl.io.json import load_json_model
import os
os.chdir(r"D:\code\github\etfl")

def extra_conc_results(model,mol_type,optim=True):
    '''extracte the proteins/mrna concentration result after model optimizing.
    enzyme.X is the concentration of the protein/mrna: mmol/gDW-1.
    PS. scaled_X is the relative concentration:fraction to the total macromolecule.
    *param:
    model:etfl formate model
    type:'enzyme' or 'mRNA'
    *return:dict_type. e.g. {'prot1':0.1,'prot2':0.2}
    e.g. enzyme_conc_result=extra_conc_results(test_model,'enzyme')
    '''
    if optim==True:
        model.optimize()
    mol_conc_result={}
    if mol_type=='enzyme':
        molecular_List=model.enzymes
    elif mol_type=='mRNA':
        molecular_List=model.mrnas
    for mol in molecular_List:
        mol_id=mol.id
        mol_conc=mol.X
        mol_conc_result[mol_id]=mol_conc
    return mol_conc_result

def gene_to_rxn(model,geneID):
    '''get the reactions from geneID
    param:
    model:etfl formate model
    geneID:the gene id
    e.g. rxnlist=gene_to_rxn(test_model,geneID)
    '''
    gene=model.genes.get_by_id(geneID)
    rxnList=list(gene.reactions)
    if len(rxnList)==0:
        print("The gene is not involved in any reaction")
        return 'None'
    elif len(rxnList)==1:
        print("The gene is involved in one reaction")
        return rxnList[0]
    else:
        print("The gene is involved in more than one reaction")
        return rxnList

# test
if __name__ == '__main__':
    # extracte the proteins concentration result after optimized
    # test_model=load_json_model("models/ME_ecoli_coreModel.json")
    yefl=load_json_model("yetfl/models/yeast8_cEFL_2584_enz_128_bins__20221115_120238.json",solver='optlang-gurobi')
    enzyme_conc_result=extra_conc_results(yefl,'enzyme')
    mrna_conc_result=extra_conc_results(yefl,'mRNA')
    # get the reaction list of the gene
    coding_geneList=list(mrna_conc_result.keys())
    for geneID in coding_geneList:
        rxnList=gene_to_rxn(yefl,geneID)


# etfl_function: 根据酶浓度得到肽链浓度
from etfl.analysis.utils import enzymes_to_peptides_conc
peps_conc=enzymes_to_peptides_conc(yefl,enzyme_conc_result)
