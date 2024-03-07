# -*- coding: utf-8 -*-
'''Load standard datasets foe strain design algorithm evaluation
'''

import pandas as pd
import os

def load_experiment_targets(product:str, data_dir='data/experiment_targets'):
    '''Load experiment targets for a specific product
    '''
    available_products = [f.replace('_exp_targets.tsv','') for f in os.listdir(data_dir) if f.endswith('_exp_targets.tsv')]
    if product not in available_products:
        print('Available products:', available_products)
        raise ValueError('The product %s is not available!' % product)
    else:
        df = pd.read_csv(os.path.join(data_dir, product+'_exp_targets.tsv'), sep='\t', index_col=0)
        return df


