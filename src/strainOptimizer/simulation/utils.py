# -*- coding: utf-8 -*-


def set_carbon_source_bounds(model, c_source, c_uptake, model_type, fixed=True):
    """Set carbon source reaction bounds based on model type.

    Args:
        fixed: True = both bounds equal c_uptake (fixed-rate simulation).
               False = one bound is 0 (open-range, e.g. yield calculation).
    """
    if model_type in ('etfl', 'GEM', 'GAN_ec'):
        if fixed:
            model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake
        else:
            model.reactions.get_by_id(c_source).bounds = -c_uptake, 0
    elif model_type == 'ecGEM':
        if fixed:
            model.reactions.get_by_id(c_source).bounds = c_uptake, c_uptake
        else:
            model.reactions.get_by_id(c_source).bounds = 0, c_uptake


def find_ngam_reaction(model):
    # 1. set of keywords to identify NGAM candidates
    keywords = ['atpm', 'maintenance', 'non-growth', 'ngam', 'atp maintenance']
    
    candidates = []
    for rxn in model.reactions:
        # Check if any of the keywords are in the reaction ID or name (case-insensitive)
        if any(key in rxn.id.lower() or key in rxn.name.lower() for key in keywords):
            candidates.append(rxn)
            
    # 2. feature validation: NGAM is typically ATP + H2O -> ADP + Pi + H+
    # Check if the reactants include ATP and the products include ADP
    final_candidates = []
    for rxn in candidates:
        mets = {met.name.lower(): coeff for met, coeff in rxn.metabolites.items()}
        # Check if it includes the key components of ATP hydrolysis
        if any('atp' in name for name in mets) and any('adp' in name for name in mets):
            final_candidates.append(rxn)

    return final_candidates