# -*- coding: utf-8 -*-

from cobra.io import read_sbml_model

def load_ecmodel(filename):
    '''load ecmodel from xml file
    parameters:
        filename: str, the path of the xml file
        return:
            model: ecmodel
            '''
    model = read_sbml_model(filename)
    return model


def ecGEM_fix_gluc(model,uptake,ex_glucID='r_1714_REV'):
    model.reactions.get_by_id(ex_glucID).bounds=uptake,uptake
    return model


