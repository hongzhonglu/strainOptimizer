# -*- coding: utf-8 -*-
# date : 2023/2/7 
# author : wangh
# file : 6.gene_modification.py
# project : etfl
'''This script is used to modify the gene expression in the model: knock-out,knock-down,overexpression.
策略：通过设置目标基因的转录/翻译反应ub/lb实现敲除/抑制，以及过表达。
'''
from etfl.io.json import load_json_model

def gene_knock_out_down_up(model,geneID,regulate_ratio):
    '''regulate the gene expression by setting the transcription/translation reaction ub/lb
    *param:
    model:etfl formate model
    geneID:the target gene id
    regulate_ratio:the ratio to change the target enzyme abundance. overexpression: >1, knock-down: 0-1, knock-out: 0
    return:modified model

    '''
    model_name=model.id
    model.optimize()
    # check weather the model has solution
    if model.solution.status!='optimal':
        print("The model has no solution")
        return 'None'
    trans_rxn=model.get_translation(geneID)
    initial_abundance=trans_rxn.flux
    rugulate_abundance = initial_abundance * regulate_ratio
    if regulate_ratio>1:
        model.logger.info("simulating %s overexpression in %s"%(geneID, model.id))
        trans_rxn.lower_bound = rugulate_abundance
    if regulate_ratio<=1:
        model.logger.info("simulating %s knock-down/out in %s"%(geneID, model_name))
        trans_rxn.upper_bound = rugulate_abundance
    # if regulate_ratio==0:
    #     model.logger.info("simulating %s knock-out in %s"%(geneID, model.id))
    #     trans_rxn.upper_bound = rugulate_abundance
    return model


if __name__ == '__main__':
    import os
    os.chdir(r'D:\code\github\etfl')
    # test_model=load_json_model("models/ME_ecoli_coreModel.json")
    yefl=load_json_model("yetfl/models/yeast8_cEFL_2584_enz_128_bins__20221228_090737.json")
    model=yefl
    sol=model.optimize()
    print("originate objective value: %s" % sol.objective_value)
    coding_geneList=[mrna.id for mrna in model.mrnas]
    for i in range(len(coding_geneList)):
        with model:
            model_change=gene_knock_out_down_up(model,geneID=coding_geneList[i],regulate_ratio=2)
            gr=model_change.slim_optimize()
            print("objective value: %s"%gr)





