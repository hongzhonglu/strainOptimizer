from strainOptimizer.io import load_model
import pandas as pd
from strainOptimizer.simulation.pprotFBA import pprotFBA_prot_conc
from strainOptimizer.strainDesign.ecFactory.ecFactory_other import genelist_to_enzymelist
from strainOptimizer.manipulation.constraint import enzyme
from strainOptimizer.analysis import optimal_yield
from etfl.optim.constraints import ModelConstraint


def change_enz_conc_bounds(model, enzID, lb=None, ub=None):
    '''change the constraint bounds in ETFL model'''
    old_constraint = model.get_constraints_of_type('ModelConstraint').get_by_id('enz_conc_'+enzID)
    expression= old_constraint.expr
    old_lb=old_constraint.constraint.lb
    old_ub=old_constraint.constraint.ub
    if lb is None:
        lb=old_lb
    if ub is None:
        ub=old_ub
    # remove the original constraint
    model.remove_constraint(old_constraint)
    # add new constraint
    model.add_constraint(kind=ModelConstraint,
                         hook=model,
                         expr=expression,
                         id_='enz_conc_'+enzID,
                         lb=lb,
                         ub=ub)
    model.repair()
    return model


# load model
model_type='etfl'
model=load_model('examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json', solver='optlang-gurobi',model_type=model_type)

c_source="r_1714"      # glucose exchange rxn
c_uptake=1
product_id='r_1589'
expYield=0.49
# load gene_enz_euva result
gene_enz_euva=pd.read_excel('examples/result/yefl_2-phenylethanol_gluc_1_ecFactory_result.xlsx',index_col=0,sheet_name='gene_enz_fva_result')
# remove rows with all NaN
gene_enz_euva=gene_enz_euva.dropna(how='all')
gene_targetList=gene_enz_euva.index.tolist()

target_enz_list, gene_enz_dict = genelist_to_enzymelist(model=model,
                                                        genelist=gene_targetList,
                                                        model_type=model_type)

# step 1. construct optimal production mutant
mutant_model = model.copy()

df_enz_bounds = pd.DataFrame(columns=['lb', 'ub'])
for geneID in gene_targetList:
    enzID = gene_enz_dict[geneID][0]
    # df_enz_bounds.loc[enzID,'lb']=target_genes_enzfva_result.loc[geneID,'prod_minprot']
    df_enz_bounds.loc[enzID, 'lb'] = 0
    df_enz_bounds.loc[enzID, 'ub'] = gene_enz_euva.loc[geneID, 'prod_max']

model.slim_optimize()
mutant_model=enzyme.ETFL_constrain_enz_conc(mutant_model,enzymes_bounds=df_enz_bounds,tol_ratio=0.01)
mutant_model.reactions.get_by_id(c_source).bounds=-1000,0
growth_id = mutant_model.growth_reaction.id

c_source_MW = 0.180156  # g/mmol
exp_gr = expYield * c_uptake * c_source_MW
# exp_gr=expYield*c_uptake
mutant_model.reactions.get_by_id(growth_id).bounds = exp_gr, exp_gr
# calculate optimal production yield,and optimal production rate
opt_prod_yield, opt_prod_rate = optimal_yield.cal_max_yield(model=mutant_model,
                                                            targetID=product_id,
                                                            c_source=c_source,
                                                            model_type=model_type)


# step 2. respectively convert each enzyme concentraction to WT-like constraints,and calculate the production yield ,and production rate
for gene in gene_targetList:
    enzID = gene_enz_dict[gene][0]
    # extract all ModelConstraint type constraints
    enz_conc_constriant=mutant_model.get_constraints_of_type('ModelConstraint').get_by_id('enz_conc_'+enzID)
    expr=enz_conc_constriant.expr
    prod_ub = enz_conc_constriant.constraint.ub
    prod_lb = enz_conc_constriant.constraint.lb
    # convert to WT-like constraint
    wt_ub = gene_enz_euva.loc[gene, 'wt_max']
    wt_lb = gene_enz_euva.loc[gene, 'wt_minprot']

    # change the enz concentration
    mutant_model=change_enz_conc_bounds(mutant_model,enzID,lb=wt_lb,ub=wt_ub)
    mod_prod_yield, mod_prod_rate = optimal_yield.cal_max_yield(model=mutant_model,
                                                                targetID=product_id,
                                                                c_source=c_source,
                                                                model_type=model_type)

    # calculate the score
    score=mod_prod_yield/opt_prod_yield+mod_prod_rate/opt_prod_rate

    print('enzyme %s, mod_prod_yield: %s, mod_prod_rate: %s, score: %s'%(enzID,mod_prod_yield,mod_prod_rate,score))

    # change back to the original constraint
    mutant_model=change_enz_conc_bounds(mutant_model,enzID,lb=prod_lb,ub=prod_ub)






