# -*- coding: utf-8 -*-
# date : 2023/3/18 
# author : wangh

def cal_max_yield(model,targetID,c_source,model_type='etfl',tol_ratio=0.001):
    '''calculate the maximum yield of the target product
    parameters:
        model: ETFL model
        targetID: str, the objective reaction ID
        c_source: str, the carbon source uptake reaction ID
    return:
        max_yield: float, the maximum yield of the target product()
    '''
    # 1. release carbon source uptake
    if model_type=='etfl':
        model.reactions.get_by_id(c_source).bounds=-1000,0
    elif model_type=='ecGEM':
        model.reactions.get_by_id(c_source).bounds=0,1000
    # 2. calculate maximum target product production
    model.objective=targetID
    model.objective_direction='max'
    max_prod=model.slim_optimize()
    if model.solver.status != 'optimal':
        return 0,0
    # 3. fix the max target product production and minimize the C source uptake
    model.reactions.get_by_id(targetID).bounds=max_prod*(1-tol_ratio),max_prod
    if model_type=='etfl':
        # model.reactions.get_by_id(c_source).bounds=-c_uptake,0
        model.objective = c_source
        model.objective_direction = 'max'
        opt_c_uptake=-model.slim_optimize()
    elif model_type=='ecGEM':
        # model.reactions.get_by_id(c_source).bounds=0,c_uptake
        model.objective=c_source
        model.objective_direction='min'
        opt_c_uptake=model.slim_optimize()

    # 4. calculate the maximum yield
    try:
        max_yield=max_prod/opt_c_uptake    # mmol/mmol carbon source
    except:
        max_yield=0
    # max_yield=abs(max_yield)

    # restore the fixed bounds
    model.reactions.get_by_id(targetID).bounds = 0,1000

    return max_yield, max_prod