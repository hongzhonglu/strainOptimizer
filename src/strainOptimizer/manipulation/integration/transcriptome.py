# -*- coding: utf-8 -*-
import pandas as pd
from .gimme import run_GIMME


def integrate_omic_data_to_ecmodel(model, omic_data, method, parameters, copy_model=False):
    '''Integrate omic data to ecmodel.
    parameters
    ----------
    model: ecGEM;
    omic_data: pd.Series
        omic data
    method: {'GIMME','soft_constraint}
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
    if method == 'GIMME':
        objective_reaction_id = parameters['objective_reaction_id']
        obj_frac = parameters['obj_frac']
        expression_threshold = parameters['expression_threshold']
        inactive_reactionID_list = run_GIMME(omic_data=omic_data, objective_reaction_id=objective_reaction_id,
                                             obj_frac=obj_frac, expression_threshold=expression_threshold)

        # remove the inactive reactions
        to_remove_rxnList = [rxn for rxn in model.reactions if rxn.id in inactive_reactionID_list]
        print(f'Removing {len(to_remove_rxnList)} reactions by {method} method')
        if copy_model:
            specific_model = model.copy()
            specific_model.remove_reactions(to_remove_rxnList)
        else:
            specific_model = model
            specific_model.remove_reactions(to_remove_rxnList)

        return specific_model
    elif method == 'soft_constraint':
        # developing...
        pass

