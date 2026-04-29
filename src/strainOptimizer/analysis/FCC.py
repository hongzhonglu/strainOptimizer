from strainOptimizer.simulation import mopa,moma

def calculate_FCC_by_abundance(protID,model,productID,c_source='r_1714_REV', c_uptake=10, growthID='r_2111',objective='r_4046',objective_direction='max',delta_conc=1):
    """
    Calculate the flux control coefficient (FCC) for a given product and growth reaction by disturb enzyme abundance.

    Args:
        model (cobra.Model): The GEM model object.
        c_source (str): The reaction ID of the carbon source uptake reaction. default is 'r_1714_REV' (glucose uptake).
        c_uptake (float): The uptake rate of the carbon source.
        productID (str): The reaction ID of the product output reaction.
        growthID (str): The reaction ID of the growth reaction.
        protID (str): The protein ID to calculate FCC for.
        objective (str): The reaction ID of the objective reaction. default is NGAM (r_4046).
        objective_direction (str): The direction of the objective reaction. 'max' or 'min'. default is 'max'.

    Returns:
        tuple: FCCg, FCCp
    """

    # calculate the reference strain by maximizing NGAM
    # ref_growth=0.2
    with model:
        model.reactions.get_by_id(c_source).bounds = (c_uptake, c_uptake)  # set uptake rate
        model.objective = growthID
        model.objective_direction='max'
        ref_growth=model.slim_optimize()/4  # set growth rate to 25% of max to allow for production
        model.objective=productID
        model.objective_direction='max'
        max_production=model.slim_optimize()
        ref_production=max_production/4
        model.reactions.get_by_id(productID).bounds = (ref_production, 1000)
        model.reactions.get_by_id(growthID).bounds = (ref_growth, 1000)  # set growth reaction bounds

        model.objective = objective  # NGAM maximize as objective
        model.objective_direction=objective_direction
        ref_solution= model.optimize()

    # calculate FCCg and FCCp
    with model:
        # overexpression for target protein
        ref_conc=ref_solution.fluxes[protID]
        new_conc=ref_conc*(1+delta_conc)  # increase protein concentration by delta_conc
        # set the protein concentration
        model.reactions.get_by_id(protID).lower_bound= new_conc

        solution=mopa(model,reference_solution=ref_solution,linear=True)
        # solution=moma(model,reference_solution=ref_solution,linear=False)

        new_growth=solution.fluxes[growthID]
        new_production=solution.fluxes[productID]

    # FCCg
    FCCg= ((new_growth-ref_growth)/ref_growth)/delta_conc
    # calculate FCCp
    FCCp=((new_production-ref_production)/ref_production)/delta_conc
    # print('growth:',new_growth,'vs',ref_growth,'production:',new_production,'vs',ref_production)

    return FCCg, FCCp


def calculate_FCC_by_kcat(protID,model,productID, c_uptake=10, growthID='r_2111',delta_kcat=1):
    '''
    Calculate the flux control coefficient (FCC) for a given product and growth reaction by disturb enzyme kcat.
    v/kcat<=protein_pool/MW
    v/(kcat*(1+delta_kcat))<=protein_pool/MW
    therefore, disturb kcat could be processed by modify the draw protein reaction coefficient
    v/kcat<=protein_pool*(1+delta_kcat)/MW
    v/kcat<=protein_pool/(MW/(1+delta_kcat))
    MW'=MW/(1+delta_kcat)

    Args:
        model (cobra.Model): The GEM model object.
        c_uptake (float): The uptake rate of the carbon source.
        productID (str): The reaction ID of the product output reaction.
        growthID (str): The reaction ID of the growth reaction.
        protID (str): The protein ID to calculate FCC for.
    '''
    c_source='r_1714_REV'  # glucose uptake reaction
    model.reactions.get_by_id(c_source).bounds = 0, c_uptake  # set uptake rate
    with model:
        model.objective=productID
        model.objective_direction='max'
        ref_production=model.slim_optimize()
        
        model.objective = growthID
        model.objective_direction='max'
        ref_growth=model.slim_optimize()
    
    # desturbe kcat
    with model:
        prot_pool=model.metabolites.get_by_id('prot_pool[c]')
        ref_mw=model.reactions.get_by_id(protID).metabolites[prot_pool]
        model.reactions.get_by_id(protID).metabolites[prot_pool]=ref_mw/(1+delta_kcat)

        model.objective = productID
        model.objective_direction='max'
        new_production=model.slim_optimize()

        model.objective = growthID
        model.objective_direction='max'
        new_growth=model.slim_optimize()

    FCCg=((new_growth-ref_growth)/ref_growth)/delta_kcat
    FCCp=((new_production-ref_production)/ref_production)/delta_kcat

    return FCCg,FCCp

    
