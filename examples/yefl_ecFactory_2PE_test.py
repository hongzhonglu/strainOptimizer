# -*- coding: utf-8 -*-
"""
2-Phenylethanol Production Strain Design using ecFactory

This script runs the ecFactory algorithm to identify genetic engineering targets
for improving 2-phenylethanol production in yeast using the ETFL model.
"""

import os
import time
import pandas as pd
import numpy as np
from pathlib import Path

from strainOptimizer.io import load_etfl_model, save_output_to_excel
from strainOptimizer.strainDesign.ecFactory import run_ecFactory
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes

# ============================================================================
# Configuration Parameters
# ============================================================================

# Model configuration
MODEL_PATH = 'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'
SOLVER = 'optlang-gurobi'

# Product and strain parameters
PRODUCT_NAME = '2-phenylethanol'
PRODUCT_ID = 'r_1589'  # 2-PE exchange reaction
C_SOURCE = 'r_1714'    # Glucose exchange reaction
C_UPTAKE = 5.0         # Glucose uptake rate (mmol/gDW/h)
GLUCOSE_MW = 0.180156  # Glucose molecular weight (g/mmol)

# Algorithm parameters
ACTION_THRESHOLDS = [0.05, 0.2, 1.1]  # Thresholds for [KO, KD, OE] classification
REMOVE_ESSENTIAL = True                # Remove essential genes from targets
STEPS = 123                            # ecFactory steps (1=FSEOF only, 12=FSEOF+EUVA, 123=full)

# Output configuration
OUTPUT_DIR = 'examples/result'
OUTPUT_FILENAME = f'yefl_{PRODUCT_NAME}_gluc_{C_UPTAKE}_ecFactory_result.xlsx'

# ============================================================================
# Model Loading and Setup
# ============================================================================

print("=" * 70)
print("Loading ETFL Model")
print("=" * 70)
model = load_etfl_model(MODEL_PATH, solver=SOLVER)
print(f"Model loaded: {MODEL_PATH}")

# Set carbon source uptake constraint
model.reactions.get_by_id(C_SOURCE).bounds = (-C_UPTAKE, -C_UPTAKE)
model.objective = model.growth_reaction.id

# Calculate maximum theoretical yield
print("\nCalculating maximum theoretical yield...")
solution = model.optimize()
if solution.status != 'optimal':
    raise RuntimeError(f"Model optimization failed with status: {solution.status}")

max_yield = solution.objective_value / (C_UPTAKE * GLUCOSE_MW)  # gDW / gGluc
print(f"Maximum theoretical yield: {max_yield:.4f} gDW/gGluc")

# Set experimental yield (49% of maximum)
exp_yield = max_yield * 0.49
alpha_limits = (0.5 * exp_yield, 2 * exp_yield)
print(f"Experimental yield: {exp_yield:.4f} gDW/gGluc")
print(f"Scanning range: [{alpha_limits[0]:.4f}, {alpha_limits[1]:.4f}] gDW/gGluc")

# ============================================================================
# Prepare Parameters for ecFactory
# ============================================================================

model_param = pd.Series({
    'targetID': PRODUCT_ID,
    'c_source': C_SOURCE,
    'c_uptake': C_UPTAKE,
    'productName': PRODUCT_NAME,
    'model_type': 'etfl'
})

# ============================================================================
# Run ecFactory Design
# ============================================================================

print("\n" + "=" * 70)
print("Running ecFactory Design Algorithm")
print("=" * 70)
print(f"Product: {PRODUCT_NAME}")
print(f"Target reaction: {PRODUCT_ID}")
print(f"Carbon source: {C_SOURCE} ({C_UPTAKE} mmol/gDW/h)")
print(f"Action thresholds: {ACTION_THRESHOLDS}")
print(f"Remove essential genes: {REMOVE_ESSENTIAL}")
print(f"Steps: {STEPS}")

start_time = time.time()
print(f"\nStart time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}")

try:
    results = run_ecFactory.run_ecFactory_design(
        model=model,
        modelParam=model_param,
        expYield=exp_yield,
        alphaLims=alpha_limits,
        action_thresholds=ACTION_THRESHOLDS,
        remove_essential=REMOVE_ESSENTIAL,
        steps=STEPS
    )
except Exception as e:
    print(f"\nError during ecFactory execution: {e}")
    raise

end_time = time.time()
elapsed_time = end_time - start_time
print(f"\nEnd time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}")
print(f"Total execution time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")

# ============================================================================
# Display Results Summary
# ============================================================================

print("\n" + "=" * 70)
print("Results Summary")
print("=" * 70)

if 'geneTable' not in results:
    raise KeyError("Results dictionary missing 'geneTable' key")

gene_table = results['geneTable']

# FSEOF results
print("\nFSEOF Results:")
oe_count = gene_table[gene_table['action'] == 'OE'].shape[0]
kd_count = gene_table[gene_table['action'] == 'KD'].shape[0]
ko_count = gene_table[gene_table['action'] == 'KO'].shape[0]
print(f"  Overexpression (OE) candidates: {oe_count}")
print(f"  Knockdown (KD) candidates: {kd_count}")
print(f"  Knockout (KO) candidates: {ko_count}")
print(f"  Total candidates: {len(gene_table)}")

# EUVA results (if available)
if 'target_priority_leval' in gene_table.columns or 'EUVA_priority_level' in gene_table.columns:
    print("\nEUVA Results:")
    priority_col = 'EUVA_priority_level' if 'EUVA_priority_level' in gene_table.columns else 'target_priority_leval'
    
    level_1 = gene_table[gene_table[priority_col] == 1].shape[0]
    level_2 = gene_table[gene_table[priority_col] == 2].shape[0]
    level_3 = gene_table[gene_table[priority_col] == 3].shape[0]
    
    print(f"  Level 1 candidates (distinct): {level_1}")
    print(f"  Level 2 candidates (overlapped): {level_2}")
    print(f"  Level 3 candidates (undistinguishable): {level_3}")

# Minimal set results (if available)
if 'minimal_candidates_set' in gene_table.columns:
    min_set_count = gene_table[gene_table['minimal_candidates_set'] == 1].shape[0]
    print(f"\nMinimal Candidate Set: {min_set_count} targets")

# ============================================================================
# Save Results to Excel
# ============================================================================

print("\n" + "=" * 70)
print("Saving Results")
print("=" * 70)

# Ensure output directory exists
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

output_path = output_dir / OUTPUT_FILENAME
print(f"Saving results to: {output_path}")

try:
    save_output_to_excel(results, str(output_path))
    print("Results saved successfully!")
except Exception as e:
    print(f"Error saving results: {e}")
    raise

print("\n" + "=" * 70)
print("Analysis Complete")
print("=" * 70)
