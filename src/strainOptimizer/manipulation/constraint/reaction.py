# -*- coding: utf-8 -*-
import math
def set_rxnbounds_1000(model):
    # set reactions with bounds more than 1000/-1000 as 1000
    for rxn in model.reactions:
        if rxn.upper_bound > 1000 or rxn.lower_bound == math.inf :
            rxn.upper_bound = 1000
        if rxn.lower_bound < -1000 or rxn.lower_bound == -math.inf:
            rxn.lower_bound = -1000

    return model