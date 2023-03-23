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