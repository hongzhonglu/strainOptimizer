# -*- coding: utf-8 -*-
# date : 2023/3/18 
# author : wangh

def cal_max_yield(model,targetID,c_source,c_uptake,model_type='etfl'):
    '''calculate the maximum yield of the target product
    parameters:
        model: ETFL model
        targetID: str, the objective reaction ID
        c_source: str, the carbon source uptake reaction ID
        c_uptake: float, the carbon source uptake rate
    return:
        max_yield: float, the maximum yield of the target product
    '''
    # 1. fix carbon source uptake
    if model_type=='etfl':
        model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake
    elif model_type=='ecGEM':
        model.reactions.get_by_id(c_source).bounds=c_uptake,c_uptake
    # 2. calculate optimal production yield
    model.objective=targetID
    model.objective_direction='max'
    max_prod=model.slim_optimize()
    if model.solver.status != 'optimal':
        return 0,0
    # 3. fix the max target product production and minimize the C source uptake
    model.reactions.get_by_id(targetID).bounds=max_prod,max_prod
    if model_type=='etfl':
        model.reactions.get_by_id(c_source).bounds=-c_uptake,0
        model.objective = c_source
        model.objective_direction = 'max'
        opt_c_uptake=-model.slim_optimize()
    elif model_type=='ecGEM':
        model.reactions.get_by_id(c_source).bounds=0,c_uptake
        model.objective=c_source
        model.objective_direction='min'
        opt_c_uptake=model.slim_optimize()

    # 4. calculate the maximum yield
    max_yield=max_prod/opt_c_uptake
    # max_yield=abs(max_yield)

    # restore the fixed bounds
    model.reactions.get_by_id(targetID).bounds = 0,1000

    return max_yield, max_prod