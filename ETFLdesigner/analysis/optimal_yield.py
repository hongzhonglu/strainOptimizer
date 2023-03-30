# -*- coding: utf-8 -*-
# date : 2023/3/18 
# author : wangh



def cal_max_yield(model,targetID,c_source,c_uptake,tol=1e-10):
    '''calculate the maximum yield of the target product

    Parameters:
        model (etfl.core.model.ETFLModel): The ETFL model
        targetID (str): The ID of the objective reaction
        c_source (str): The ID of the carbon source uptake reaction
        c_uptake (float): The carbon source uptake rate
        tol (float, optional): The tolerance for the optimization. Defaults to 1e-10.

    Returns:
        max_yield (float): The maximum yield of the target product
        max_prod (float): The maximum production of the target product
    '''
    # 1. fix carbon source uptake
    model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake
    # 2. calculate optimal production yield
    model.objective=model.reactions.get_by_id(targetID)
    model.objective_direction='max'
    max_prod=model.slim_optimize()
    # 3. fix the max target product production and minimize the C source uptake
    model.reactions.get_by_id(targetID).bounds=max_prod-tol,max_prod+tol
    model.objective=c_source
    model.objective_direction='min'
    opt_c_uptake=-model.slim_optimize()
    max_yield=max_prod/opt_c_uptake

    return max_yield, max_prod
