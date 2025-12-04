# -*- coding: utf-8 -*-
'''Run Thermodynamics-based Flux Analysis(TFA).
reference:https://github.com/EPFL-LCSB/pytfa'''
from pytfa.io import load_thermoDB,read_lexicon,read_compartment_data,annotate_from_lexicon,apply_compartment_data
import pandas as pd
import pytfa
from strainOptimizer.manipulation.constraint.reaction import set_rxnbounds_1000

def curate_lexicon(lexicon):
    ix = pd.Series(lexicon.index)
    ix = ix.apply(lambda s: str.replace(s, '-', '__'))
    ix = ix.apply(lambda s: '_' + s if s[0].isdigit() else s)
    lexicon.index = ix
    return lexicon


def build_tfa_model(model,objective='r_2111',model_type='ecGEM',relax_threshold=0.5,yeastGEM=True,compartment_data_path=None):
    '''
    Build a TFA model from a GEM/ecGEM model.
    parameters:
        model: GEM/ecGEM model
        objective: target reaction ID
        model_type: 'GEM' or 'ecGEM'
        yeastGEM: Check if the model is built from yeastGEM(https://github.com/SysBioChalmers/yeast-GEM), default True,
                    if False, the model will be built the tfa model by eQuilibrator in pyTFA.
        relax_threshold: relaxation threshold for the TFA model, default 0.5. It means if the tfa optimized value if less
                        than 0.5*fba_value, the relax process will be triggered.
        compartment_data_path: the path of the compartment data file, default None. If None, the S. cerevisiae specific
                        compartment data will be loaded from the default path.
    return:
        tfa_model: pytfa.ThermoModel
    '''
    if model_type not in ['GEM','ecGEM']:
        raise ValueError('model_type should be either "GEM" or "ecGEM"')

    model=set_rxnbounds_1000(model)

    # build the tfa model for yeastGEM-based model
    if yeastGEM is True:
        if model_type=='GEM':
            if 'v9' in model.id:
                lexicon_path=r'data/yeast_thermo/yeast9_lexicon.csv'
            else:
                lexicon_path=r'data/yeast_thermo/yeast_lexicon.csv'
        elif model_type=='ecGEM':
            lexicon_path=r'data/yeast_thermo/ecYeast_lexicon.csv'
        thermoDB_path=r'data/yeast_thermo/thermo_data.thermodb'
        compartment_data_path=r'data/yeast_thermo/yeast_compartment_data.json'

        # load data
        thermo_data=load_thermoDB(path=thermoDB_path)
        lexicon=read_lexicon(lexicon_path)
        lexicon=curate_lexicon(lexicon)
        compartment_data=read_compartment_data(compartment_data_path)

        # initialize the tfa model
        model.objective=objective
        tfa_model=pytfa.ThermoModel(thermo_data,model)
        tfa_model.name=model.name

        # add thermo information
        annotate_from_lexicon(tfa_model, lexicon)
        apply_compartment_data(tfa_model, compartment_data)

        # prepare and convert the model
        tfa_model.prepare()
        tfa_model.convert()

    else:
        # build the tfa model by eQuilibrator in pyTFA
        from pytfa.thermo.equilibrator import build_thermo_from_equilibrator
        thermo_data = build_thermo_from_equilibrator(model)
        if compartment_data_path is None:
            compartment_data_path = r'data/yeast_thermo/yeast_compartment_data.json'
        compartment_data = read_compartment_data(compartment_data_path)

        tfa_model = pytfa.ThermoModel(thermo_data, model)
        apply_compartment_data(compartment_data=compartment_data, tmodel=tfa_model)
        tfa_model.prepare()
        tfa_model.convert()

    # check wheather the model is over-constrained by thermodynamics, if so, run a relax procedure
    fba_value=model.slim_optimize()
    tfa_value=tfa_model.optimize().objective_value

    if tfa_value<fba_value*relax_threshold:
        from pytfa.optim.relaxation import relax_dgo

        tfa_model.reactions.get_by_id(objective).lower_bound = 0.5 * fba_value
        relaxed_model, slack_model, relax_table = relax_dgo(tfa_model)
        print('The model is over-constrained by thermodynamics, run a relax procedure.')
        print(f'{len(relax_table)} reactions are relaxed:')
        print(relax_table)

        return relaxed_model
    else:
        return tfa_model


def tfa(model, targetID, c_source, c_uptake, model_type='ecGEM', direction='max',yeastGEM=True):
    if model_type not in ['etfl', 'ecGEM','GEM']:
        raise ValueError('model_type should be either "etfl" or "ecGEM" or "GEM"')

    # set the carbon source
    if c_uptake is None:
        c_uptake = 1
    if model_type == 'etfl':
        model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake
    elif model_type == 'ecGEM':
        model.reactions.get_by_id(c_source).bounds = c_uptake, c_uptake
    elif model_type == 'GEM':
        model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake

    # if the model GEM or ecGEM, build a TFA model
    if model_type in ['GEM','ecGEM']:
        model=build_tfa_model(model,objective='r_2111',model_type=model_type,yeastGEM=yeastGEM)

    # run TFA
    model.objective = targetID
    model.objective_direction = direction
    sol = model.optimize()

    return sol
