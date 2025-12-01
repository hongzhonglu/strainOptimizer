# -*- coding: utf-8 -*-
"""
Comprehensive Example for Strain Design Workflow Engine

This script demonstrates various usage patterns for the strainOptimizer workflow engine.
Based on real usage patterns from 1.ecFactory_ecYeast_test.py, this shows how to:
1. Run ecFactory with specialized ecGEM models
2. Configure parameters using plain dictionaries
3. Update parameters after instantiating the engine
"""

from strainOptimizer.strainDesign.workflow_engine import (
    strainOptimizer_engine,
    WorkflowParameters,
)


def example_real_artemisinic_acid_workflow():
    """
    Example 1: Real workflow for artemisinic acid production using new framework
    Based on the actual ecFactory_ecYeast_test.py workflow
    Demonstrates the dictionary-based parameter structure:
    - model: Model-related parameters
    - strain: Strain phenotype parameters
    - algorithm: Algorithm control parameters
    """
    print("=== Example 1: Real Artemisinic Acid Production Workflow ===")

    # Method 1: Using dictionaries (recommended)
    # Model parameters - model file and solver settings
    model_params = {
        'model_path': 'examples/models/yeast/ecGEM_atemisinic.xml',
        'model_type': 'ecGEM',
        'solver': 'optlang-gurobi',
        'growth_id': 'r_2111',
        # 'total_enzymes': 0.1
    }

    # Strain parameters - target product and growth conditions
    strain_params = {
        'target_id': 'r_1589',
        'product_name': '2-phenylethanol',
        'c_source': 'r_1714_REV',  # glucose exchange reaction
        'c_uptake': 1.0,  # glucose uptake rate (mmol/gDW/h)
    }

    # Algorithm control parameters - workflow and output settings
    algorithm_params = {
        'design_algorithm': 'ecFactory',
        'remove_essential': False,
        'output_directory': './results',
        'save_results': True,
        # 'only_final_result': True,
        # Note: ecFactory-specific parameters like steps, action_thresholds, etc.
        # would need to be added to AlgorithmControl if they're used
    }

    # Create WorkflowParameters using the three-level structure
    params = WorkflowParameters(
        model=model_params,
        strain=strain_params,
        algorithm=algorithm_params
    )

    try:
        # Create workflow engine using the new framework
        engine = strainOptimizer_engine(params)
        
        print(f"Engine created for {params.strain['product_name']} production")
        print(f"Target reaction: {params.strain['target_id']}")
        print(f"Carbon source: {params.strain['c_source']}")
        print(f"Model type: {params.model['model_type']}")
        print(f"Algorithm: {params.algorithm['design_algorithm']}")
        
        # Load model
        model=engine.load_model()
        
        # Get model information
        model_info = engine.get_model_info()
        print(f"\nModel info: {model_info}")
        
        # Run the design workflow
        print("\nRunning strain design workflow...")
        results = engine.run_design()
        
        # Analyze results exactly like the original workflow
        if results and 'geneTable' in results:
            gene_table = results['geneTable']
            print(f"\nFSEOF results:")
            print(f"  OE: {len(gene_table[gene_table['action'] == 'OE'])}")
            print(f"  KD: {len(gene_table[gene_table['action'] == 'KD'])}")  
            print(f"  KO: {len(gene_table[gene_table['action'] == 'KO'])}")
            
            # Check EUVA results like original
            if 'target_priority_level' in gene_table.columns:
                print(f"\nEUVA results:")
                for level in [1, 2, 3]:
                    count = len(gene_table[gene_table['target_priority_level'] == level])
                    print(f"  Level {level} candidates: {count}")
        
        # Get summary
        summary = engine.get_results_summary()
        print(f"\nDesign workflow completed successfully!")
        print(f"Summary: {summary}")
        
        return engine, results
        
    except FileNotFoundError as e:
        print(f"Model file not found: {e}")
        print("Note: This uses example paths - update with real model files for actual usage")
        return None, None
    except Exception as e:
        print(f"Workflow failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None


def example_quick_dict_setup():
    """
    Example 2: Creating parameters directly with dictionaries for a different scenario
    """
    print("\n=== Example 2: Quick Dict-Based Setup ===")
    
    params = WorkflowParameters(
        model={
            'model_path': 'examples/models/yeast/ecGEM_atemisinic.xml',
            'model_type': 'ecGEM',
        },
        strain={
            'target_id': 'r_2111',
            'product_name': 'Artemisinic acid',
            'c_uptake': 0.8,
        },
        algorithm={
            'design_algorithm': 'iBridge',
            'save_results': True,
        }
    )
    
    print("Created parameter set with:")
    print(f"  Model: {params.model['model_type']} @ {params.model['model_path']}")
    print(f"  Strain target: {params.strain['target_id']} ({params.strain['product_name']})")
    print(f"  Algorithm: {params.algorithm['design_algorithm']} (save_results={params.algorithm['save_results']})")
    
    return params


def example_parameter_updates():
    """
    Example 3: Updating parameters after creation
    """
    print("\n=== Example 3: Updating Parameters ===")
    
    # Create initial parameters
    params = WorkflowParameters(
        model={
            'model_path': 'examples/models/yeast/ecGEM_atemisinic.xml',
            'model_type': 'ecGEM',
        },
        strain={
            'target_id': 'r_1589',
            'product_name': '2-phenylethanol',
        }
    )
    
    print(f"Initial algorithm: {params.algorithm['design_algorithm']}")
    print(f"Initial save_results: {params.algorithm['save_results']}")
    
    # Update parameters using the update_parameters method
    params.update_parameters(
        design_algorithm='iBridge',
        save_results=True,
        c_uptake=2.0
    )
    
    print(f"Updated algorithm: {params.algorithm['design_algorithm']}")
    print(f"Updated save_results: {params.algorithm['save_results']}")
    print(f"Updated c_uptake: {params.strain['c_uptake']}")
    
    return params


def main():
    """
    Main function demonstrating real workflow patterns based on ecFactory_ecYeast_test.py
    Shows different ways to use the dictionary-based parameter structure
    """
    # Example 1: Using dictionaries (recommended for simplicity)
    example_real_artemisinic_acid_workflow()
    
    # Example 2: Creating parameters directly with dictionaries
    example_quick_dict_setup()
    
    # Example 3: Updating parameters
    example_parameter_updates()
    
    print("\n" + "="*60)
    print("All examples completed!")
    print("="*60)

if __name__ == "__main__":
    main()
