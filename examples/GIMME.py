# -*- coding: utf-8 -*-
import pandas as pd
import cobra
import re
from strainOptimizer.io import load_model

from troppo.methods_wrappers import ModelBasedWrapper
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties


def gene_to_rxn_calculator(gene_data,rxn,ignore_none_value=True):
    '''calculate the related transcriptional level for each reaction according to the GPR rules.
    calculating
    rules: max for OR, min for AND.
    parameters
    ----------
    gene_data: pd.Series
    rxn: cobra.Reaction
        reaction object
    '''
    # get related genes
    geneIDlist=[gene.id for gene in rxn.genes]
    # only keep the genes that are in the gene_data
    geneIDlist=[gene for gene in geneIDlist if gene in gene_data.index]

    # if there is no gene in the gene_data, return 1000000
    if len(geneIDlist)==0:
        if ignore_none_value:
            return 1000000
        else:
            raise ValueError(f"All genes in the reaction {rxn.id} are not found in the gene_data")

    rxn_gene_dataMatrix=gene_data[geneIDlist]

    # calculate the related copy number for each reaction according to the GPR rules.
    gpr=rxn.gene_reaction_rule

    # if there is no and in gpr, then the copy number is the max of related genes
    if 'and' not in gpr:
        # rxn_data=rxn_gene_dataMatrix.max(axis=0) # max for OR
        rxn_value=rxn_gene_dataMatrix.sum(axis=0)   # sum for OR
    # if there is no or in gpr, then the copy number is the min of related genes
    elif 'or' not in gpr:
        rxn_value=rxn_gene_dataMatrix.min(axis=0)
    # if there is both and and or in gpr, then firstly, divede the gpr into several parts according to or,
    # then calculate the copy number for each part by min, finally, the copy number is the max of these parts.
    else:
        df_complexes_dataMatrix=pd.DataFrame()
        complexlist=gpr.split('or')
        complexes_values=[]
        for complex in complexlist:
            # remove the space in both sides
            complex=complex.strip()
            # remove the bracket in both sides
            complex=complex.strip('(')
            complex=complex.strip(')')
            geneIDlist=complex.split('and')
            # remove the space in both sides
            geneIDlist=[geneID.strip() for geneID in geneIDlist]
            complex_data=rxn_gene_dataMatrix[geneIDlist].min(axis=0)
            complexes_values.append(complex_data)
        # sum the complexes data
        rxn_value=pd.Series(complexes_values).sum()
    return rxn_value


def gene_to_rxn_scores(gene_data,model):
    rxn_scores={}
    for rxn in model.reactions:
        rxn_scores[rxn.id]=gene_to_rxn_calculator(gene_data=gene_data,rxn=rxn)
    return rxn_scores


def run_GIMME(omic_data,objective_reaction_id,gem_model_file_path=r'data/yeast-GEM.xml',obj_frac=0.8,expression_threshold=12):
    '''Run GIMME algorithm to integrate transcriptome data with GEM model.
    parameters
    ----------
    omic_data: pd.Series
        transcriptome data
    gem_model_file_path: str
        path to the classical GEM model
    objective_reaction_id: str
        id of the objective reaction
    obj_frac: float
        fraction of the objective reaction
    expression_threshold: float
        expression threshold
    '''
    # ignore the UserWarning in cobamp package
    import warnings
    warnings.filterwarnings("ignore", category=UserWarning, module="cobamp")

    # Load the GEM model (replace with your model path)
    gem_model = cobra.io.read_sbml_model(gem_model_file_path)
    # wrap the model as troppo model
    model_wrapper = ModelBasedWrapper(model=gem_model, ttg_ratio=9999)

    # calculate the scores for each reaction
    scores = gene_to_rxn_scores(omic_data,gem_model)

    # Get the index of the biomass reaction in the model. This will be used as objective for the GIMME algorithm.
    idx_objective = model_wrapper.model_reader.r_ids.index(objective_reaction_id)

    # Create the properties for the GIMME algorithm.
    expression_threshold = float(expression_threshold)
    reactionID_list=model_wrapper.model_reader.r_ids
    metaboliteID_list=model_wrapper.model_reader.m_ids
    properties = GIMMEProperties(exp_vector=[v for k, v in scores.items()], obj_frac=obj_frac,
                                 objectives=[{idx_objective: 1}],
                                 preprocess=True,flux_threshold=expression_threshold,solver='GUROBI',
                                 reaction_ids=reactionID_list,
                                 metabolite_ids=metaboliteID_list)

    # Run the GIMME algorithm.
    gimme = GIMME(S=model_wrapper.S, lb=model_wrapper.lb, ub=model_wrapper.ub, properties=properties)

    active_reactionIndex_list = gimme.run()
    active_reactionID_list=[reactionID_list[i] for i in active_reactionIndex_list]

    inactive_reactionID_list=list(set(reactionID_list)-set(active_reactionID_list))

    return inactive_reactionID_list


def integrate_omic_data_to_ecmodel(model,omic_data,method,parameters,copy_model=False):
    '''Integrate omic data to ecmodel.
    parameters
    ----------
    model: ecGEM;
    omic_data: pd.Series
        omic data
    method: {'GIMME'}
        method to integrate omic data
    parameters: dict
        parameters for the integration method
        GIMME:
            objective_reaction_id: str
                id of the objective reaction
            obj_frac: float
                fraction of the objective reaction
            expression_threshold: float
                expression threshold
    '''
    if method=='GIMME':
        objective_reaction_id=parameters['objective_reaction_id']
        obj_frac=parameters['obj_frac']
        expression_threshold=parameters['expression_threshold']
        inactive_reactionID_list=run_GIMME(omic_data=omic_data,objective_reaction_id=objective_reaction_id,
                                         obj_frac=obj_frac,expression_threshold=expression_threshold)
        
        # remove the inactive reactions
        to_remove_rxnList=[rxn for rxn in model.reactions if rxn.id in inactive_reactionID_list]
        print(f'Removing {len(to_remove_rxnList)} reactions by {method} method')
        if copy_model:
            specific_model=model.copy()
            specific_model.remove_reactions(to_remove_rxnList)
        else:
            specific_model=model
            specific_model.remove_reactions(to_remove_rxnList)

        return specific_model
    elif method=='soft_constraint':
        # developing...
        pass


if __name__ == "__main__":
    # load the ecmodel
    from strainOptimizer.io import load_model
    model = load_model(r'examples/models/yeast/ecYeastGEM_batch.xml',model_type='ecGEM')

    # load the omic data
    df_transcriptomics = pd.read_csv(r'sclareol_cell_factory_design\data\sce_sclareol_gene_fpkm.tsv', index_col=0,
                                     sep='\t')
    samples = ['SCX42-1', 'SCX42-2', 'SCX42-3']
    df_scx42 = df_transcriptomics[samples].mean(axis=1)

    # parameters for GIMME
    parameters = {
        'objective_reaction_id': 'r_2111',
        'obj_frac': 0.8,
        'expression_threshold': 12
    }

    # run the integration
    specific_model = integrate_omic_data_to_ecmodel(model=model, omic_data=df_scx42, method='GIMME', parameters=parameters,copy_model=True)
    specific_model.slim_optimize()
