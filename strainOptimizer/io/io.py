from ETFL import load_etfl_model
from ecModel import load_ecmodel

def load_model(filename:str , model_type: ['etfl','ecGEM'] ,solver=None):
    '''load ecModel/ETFL model
    parameters:
        filename: str, the path of the model file
        model_type: str, 'etfl' or 'ecGEM'
        solver: str, the solver used to solve the model
        return:
            model: ecModel/ETFL model
            '''
    if model_type=='ec':
        model=load_ecmodel(filename)
    elif model_type=='etfl':
        model=load_etfl_model(filename,solver=solver)
    return model
