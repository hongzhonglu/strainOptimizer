# -*- coding: utf-8 -*-
from etfl.io.json import load_json_model

def load_etfl_model(filename,solver=None):
    '''load ETFL model from json file
    parameters:
        filename: str, the path of the json file
        return:
            model: ETFL model
            '''
    model = load_json_model(filename,solver=solver)
    return model