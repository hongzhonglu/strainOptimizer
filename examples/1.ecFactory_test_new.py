# -*- coding: utf-8 -*-
"""
Example script for running ecFactory strain design using the new workflow engine framework.

This script demonstrates how to use the strainOptimizer_engine to run ecFactory design
for artemisinic acid production in ecYeastGEM model.
"""

import pandas as pd
import numpy as np
import time
import cobra
from strainOptimizer.strainDesign.workflow_engine import (
    strainOptimizer_engine,
    WorkflowParameters,
)
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes

# Set tolerance
cobra.Configuration().tolerance = 1e-9

# Product parameters
productParam_dict = {
    # Add more products as needed
    '2-phenylethanol': {
        'product_name': '2-phenylethanol',
        'target_id': 'r_1589',
        'model_filepath': 'examples/models/yeast/ecYeastGEM_batch.xml'
    },
}

# Model parameters
modelParams_dict = {
    'ecGEM': {
        'model_type': 'ecGEM',
        'c_source': 'r_1714_REV',  # glucose exchange rxn
        'c_uptake': 1,
        'growth_id': 'r_2111',
    },


}

# Process each product and model combination
for model_key, modelParam in modelParams_dict.items():
    if model_key == 'etfl':
        continue
    
    model_consistency_result = pd.DataFrame(index=['l1', 'l2', 'l3'])
    
    for product_key, productParam in productParam_dict.items():
        print(f'Product: {product_key}, Model: {model_key}')
        
        # Merge product and model parameters
        modelParams = {**modelParams_dict[model_key], **productParam}
        
        # ============================================================
        # Step 1: Create WorkflowParameters
        # ============================================================
        parameters = WorkflowParameters(
            model={
                'model_path': modelParams['model_filepath'],
                'model_type': modelParams['model_type'],
                'solver': 'optlang-gurobi',
                'growth_id': modelParams['growth_id'],
            },
            strain={
                'target_id': modelParams['target_id'],
                'product_name': modelParams['product_name'],
                'c_source': modelParams['c_source'],
                'c_uptake': modelParams['c_uptake'],

                'substrate_MW': 0.180156,  # g/mmol, glucose molecular weight
            },
            algorithm={
                'design_algorithm': 'ecFactory',
                'simulation_method': 'ppfba',
                'remove_essential': True,
                'action_thresholds': [0.05, 0.3, 1.1],  # KO, KD, OE thresholds
                'steps': 123,
                'save_results': False,
                'output_directory': './results',
            }
        )
        
        # ============================================================
        # Step 2: Initialize workflow engine
        # ============================================================
        engine = strainOptimizer_engine(parameters)
        
        # ============================================================
        # Step 3: Load model
        # ============================================================
        print('\nLoading model...')
        model = engine.load_model()
        
        
        # ============================================================
        # Step 7: Run strain design
        # ============================================================
        
        start_time = time.time()
        print(f'Start time: {start_time}')
        
        # Run design
        final_result = engine.run_design()
        
        end_time = time.time()
        print(f'End time: {end_time}')
        print(f'Time cost: {end_time - start_time:.2f} seconds')
        