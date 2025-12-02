from .etfl.io.json import load_json_model
from cobra.io import read_sbml_model
import pandas as pd
import numpy as np

def load_etfl_model(filename,solver=None):
    '''load ETFL model from json file
    parameters:
        filename: str, the path of the json file
        return:
            model: ETFL model
            '''
    model = load_json_model(filename,solver=solver)
    return model



def load_ecmodel(filename):
    '''load ecmodel from xml file
    parameters:
        filename: str, the path of the xml file
        return:
            model: ecmodel
            '''
    # model = read_sbml_ec_model(filename)
    model=read_sbml_model(filename)
    return model


def load_model(filename:str , model_type: ['etfl','ecGEM','gem'] ,solver=None):
    '''load ecModel/ETFL model
    parameters:
        filename: str, the path of the model file
        model_type: str, 'etfl' or 'ecGEM'
        solver: str, the solver used to solve the model
        return:
            model: ecModel/ETFL model
            '''
    if model_type=='ecGEM' or model_type=='GAN_ec':
        model=load_ecmodel(filename)
    elif model_type=='etfl':
        model=load_etfl_model(filename,solver=solver)
    elif model_type=='gem':
        model=read_sbml_model(filename)
    # set tolerance as 1e-9
    model.tolerance = 1e-9
    return model



def save_output_to_excel(result,file_path:str):
    '''save the strainOptimizer output to excel file
    parameters:
        result: dict, output of strainOptimizer
        file_path: str, the path of the output file

            '''
    # save results into excel file
    for key in result.keys():
        if isinstance(result[key], np.ndarray):
            result[key] = pd.Series(result[key])
        # 检查是否为dict
        if isinstance(result[key], dict):
            result[key] = pd.Series(result[key])
    with pd.ExcelWriter(file_path) as writer:
        for key in result.keys():
            result[key].to_excel(writer, sheet_name=key)
    return None