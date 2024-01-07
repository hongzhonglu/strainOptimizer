from etfl.io.json import load_json_model
from cobra.io import read_sbml_model

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


def load_model(filename:str , model_type: ['etfl','ecGEM'] ,solver=None):
    '''load ecModel/ETFL model
    parameters:
        filename: str, the path of the model file
        model_type: str, 'etfl' or 'ecGEM'
        solver: str, the solver used to solve the model
        return:
            model: ecModel/ETFL model
            '''
    if model_type=='ecGEM':
        model=load_ecmodel(filename)
    elif model_type=='etfl':
        model=load_etfl_model(filename,solver=solver)
    return model
