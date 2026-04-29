# -*- coding: utf-8 -*-
"""
Strain Design Workflow Engine

This module provides a comprehensive workflow engine for strain design operations,
integrating various strain design algorithms and providing a unified interface.
"""

import copy
import pandas as pd
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
from pathlib import Path

from ..io import load_model, save_output_to_excel
from .ecFactory import run_ecFactory_design
from .ecFactory.ecfseof import run_ecFSEOF
from .iBridge import run_iBridge_design
from ..analysis.dataset import gene_id_to_name

DEFAULT_MODEL_PARAMETERS: Dict[str, Any] = {
    'model_path': None,
    'model_type': None,  # 'etfl', 'ecGEM', 'gem', 'GAN_ec'
    'solver': 'optlang-gurobi',
    'growth_id': None,
}

DEFAULT_STRAIN_PARAMETERS: Dict[str, Any] = {
    'target_id': None,
    'product_name': 'PRODUCT',
    'c_source': 'r_1714_REV',  # Default glucose exchange reaction of ecYeastGEM
    'c_uptake': 1.0,  # Default glucose uptake rate 1 mmol/gDW/h
    'substrate_MW': 0.180156,  # g/mmol, glucose molecular weight
    'expYield': None,
    'transcriptome': None
}

DEFAULT_ALGORITHM_PARAMETERS: Dict[str, Any] = {
    'design_algorithm': 'ecFactory',  # 'ecFactory', 'iBridge', 'optForce'
    'simulation_method': 'ppfba',  # 'ppfba', 'pfba', 'moma', or 'mopa'
    'remove_essential': False,
    'output_directory': './results',
    'save_results': False,
    'experimental_yield': None,
    'scanning_range': None,
    'action_thresholds': [0.05, 0.3, 1.05],
    'steps': 123,
    'calculate_fcc': False,  # whether to calculate FCC scores for final result genes
}


@dataclass
class WorkflowParameters:
    """Parameters for strain design workflow with three dictionaries"""
    model: Dict[str, Any] = field(default_factory=dict)
    strain: Dict[str, Any] = field(default_factory=dict)
    algorithm: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Ensure dictionaries are initialized with defaults"""
        provided_model = self.model or {}
        provided_strain = self.strain or {}
        provided_algorithm = self.algorithm or {}

        self.model = {**DEFAULT_MODEL_PARAMETERS, **provided_model}
        self.strain = {**DEFAULT_STRAIN_PARAMETERS, **provided_strain}
        self.algorithm = copy.deepcopy(DEFAULT_ALGORITHM_PARAMETERS)
        self.algorithm.update(provided_algorithm)
    
    def update_parameters(self, **kwargs) -> None:
        """
        Update workflow parameters. Supports passing section dicts directly
        (model/strain/algorithm) or individual keys that match existing entries.
        """
        for key, value in kwargs.items():
            if key in ('model', 'strain', 'algorithm') and isinstance(value, dict):
                getattr(self, key).update(value)
                continue

            if key in self.model:
                self.model[key] = value
            elif key in self.strain:
                self.strain[key] = value
            elif key in self.algorithm:
                self.algorithm[key] = value
            else:
                # default to algorithm dict for unknown keys
                self.algorithm[key] = value
    
    def get_parameters(self) -> Dict[str, Any]:
        """Get workflow parameters as a flat dictionary"""
        params = {}
        params.update({f"model.{k}": v for k, v in self.model.items()})
        params.update({f"strain.{k}": v for k, v in self.strain.items()})
        params.update({f"algorithm.{k}": v for k, v in self.algorithm.items()})
        return params


class strainOptimizer_engine:
    """
    A comprehensive workflow engine for strain design operations.
    
    This class provides a unified interface for running various strain design
    algorithms and managing the entire workflow from model loading to result analysis.
    """
    
    def __init__(self, parameters: WorkflowParameters):
        """
        Initialize the workflow engine with parameters.
        
        Args:
            parameters: WorkflowParameters object containing all workflow settings
        """
        self.parameters = parameters
        self.model = None
        self.all_results = {}
        self.final_result = None
        
        # Validate parameters
        self._validate_parameters()
    
    def load_model(self) -> None:
        """Load the metabolic model"""
        self._load_model()
        return self.model
        
    
    def _validate_parameters(self) -> None:
        """Validate workflow parameters"""

        # Validate required model parameters
        if not self.parameters.model.get('model_path'):
            raise ValueError("Required parameter 'model.model_path' is missing or None")
        if not self.parameters.model.get('model_type'):
            raise ValueError("Required parameter 'model.model_type' is missing or None")
        # Validate model type
        valid_model_types = ['etfl', 'ecGEM', 'gem', 'GAN_ec']
        if self.parameters.model.get('model_type') not in valid_model_types:
            raise ValueError(f"Invalid model_type. Must be one of: {valid_model_types}")
        # Validate model file existence
        model_path = self.parameters.model.get('model_path')
        if not Path(model_path).exists():
            raise FileNotFoundError(f"Model file not found: {model_path}")

        # Validate required strain parameters
        if not self.parameters.strain.get('target_id'):
            raise ValueError("Required parameter 'strain.target_id' is missing or None")

        # Validate required algorithm parameters
        valid_algorithms = ['ecFactory', 'iBridge', 'ecFSEOF']
        if self.parameters.algorithm.get('design_algorithm') not in valid_algorithms:
            raise ValueError(f"Invalid design_algorithm. Must be one of: {valid_algorithms}")
        
        print("Parameters validated successfully")
    
    def _load_model(self) -> None:
        """Load the metabolic model"""
        try:
            model_type = self.parameters.model.get('model_type')
            model_path = self.parameters.model.get('model_path')
            solver = self.parameters.model.get('solver')
            print(f"Loading {model_type} model from {model_path}")
            self.model = load_model(
                filename=model_path,
                model_type=model_type,
                solver=solver
            )
        except Exception as e:
            print(f"Failed to load model: {str(e)}")
            raise
    
    def run_design(self) -> Dict[str, Any]:
        """
        Run the strain design workflow.
        
        Returns:
            Dictionary containing design results
        """
        # check model
        if self.model is None:
            raise ValueError("Model not loaded")
        design_algorithm = self.parameters.algorithm.get('design_algorithm')
        print(f"Starting strain design workflow using {design_algorithm}")
        
        try:
            if design_algorithm == 'ecFactory':
                results = self._run_ecFactory()
            elif design_algorithm == 'iBridge':
                results = self._run_iBridge()
            elif design_algorithm == 'ecFSEOF':
                results = self._run_ecFSEOF()
            else:
                raise ValueError(f"Unsupported design algorithm: {design_algorithm}")
            
            self.all_results.update(results)
            # add gene names — only available for S. cerevisiae s288c strain models
            try:
                results['geneTable'].loc[:, 'gene_name'] = gene_id_to_name(results['geneTable'].index)
                results['level 1 result'].loc[:, 'gene_name'] = gene_id_to_name(results['level 1 result'].index)
                if 'level 2 result' in results:
                    results['level 2 result'].loc[:, 'gene_name'] = gene_id_to_name(results['level 2 result'].index)
                if 'level 3 result' in results:
                    results['level 3 result'].loc[:, 'gene_name'] = gene_id_to_name(results['level 3 result'].index)
            except Exception:
                pass
            self.final_result = results['geneTable']

            # Build FCC-filtered result when FCC scores are available
            if 'fcc_result' in results:
                results['fcc_filtered_result'] = self._fcc_filter(results['geneTable'])
                self.all_results['fcc_filtered_result'] = results['fcc_filtered_result']
                print(f"FCC filtered result: {len(results['fcc_filtered_result'])} genes retained")

            
            # Save results if requested
            if self.parameters.algorithm.get('save_results'):
                self._save_results()
            
            print("Strain design workflow completed successfully")
            return self.final_result
            
        except Exception as e:
            print(f"Workflow failed: {str(e)}")
            raise
    
    def _run_ecFactory(self) -> Dict[str, Any]:
        """Run ecFactory design algorithm"""
        # 1. check and prepare essential parameters for ecFactory workflow: expYield
        if self.parameters.algorithm.get('design_algorithm') != 'ecFactory':
            raise ValueError("Design algorithm must be 'ecFactory'")

        # 2. run ecFactory
        print("Running ecFactory algorithm")

        # Note: These parameters may need to be added to AlgorithmControl if they're algorithm-specific
        results = run_ecFactory_design(
            model=self.model,
            parameters=self.parameters
        )
        
        return results
    
    def _run_iBridge(self) -> Dict[str, Any]:
        """Run iBridge design algorithm"""
        print("Running iBridge algorithm")
        
        # Prepare parameters for iBridge
        # Note: This is a placeholder - actual implementation depends on iBridge interface
        results = run_iBridge_design(
            model=self.model,
            target_reaction=self.parameters.strain['target_id'],
            # Add other iBridge-specific parameters as needed
        )
        
        return results
    
    def _run_ecFSEOF(self) -> Dict[str, Any]:
        """Run ecFSEOF design algorithm"""
        print("Running ecFSEOF algorithm")
        
        # Prepare parameters for ecFSEOF
        results = run_ecFSEOF(
            model=self.model,
            parameters=self.parameters
        )
        
        return results
    
    def _fcc_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        """Return rows consistent with FCC scores: OE→FCCp>0, KD/KO→FCCp<0."""
        if 'FCCp' not in df.columns:
            return pd.DataFrame(columns=df.columns)
        df = df.copy()
        df['FCCp'] = df['FCCp'].apply(lambda x: 0 if pd.notna(x) and abs(x) < 1e-9 else x)
        oe_mask    = (df['action'] == 'OE')  & (df['FCCp'] > 0)
        kd_ko_mask = (df['action'].isin(['KD', 'KO'])) & (df['FCCp'] < 0)
        return df[oe_mask | kd_ko_mask]

    def _save_results(self) -> None:
        """Save workflow results to files"""
        if self.final_result is None:
            print("No results to save")
            return

        # Create output directory if it doesn't exist
        output_dir = Path(self.parameters.algorithm.get('output_directory', './results'))
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate output filename prefix
        filename_prefix = f"{self.parameters.strain.get('product_name')}_design_results_{self.parameters.algorithm.get('design_algorithm')}"

        try:
            self.final_result.to_excel(output_dir / (filename_prefix + '_combined_result.xlsx'))
            self.all_results['level 1 result'].to_excel(output_dir / (filename_prefix + '_level_1_result.xlsx'))
            if 'level 2 result' in self.all_results:
                self.all_results['level 2 result'].to_excel(output_dir / (filename_prefix + '_level_2_result.xlsx'))
            if 'level 3 result' in self.all_results:
                self.all_results['level 3 result'].to_excel(output_dir / (filename_prefix + '_level_3_result.xlsx'))
            if 'fcc_filtered_result' in self.all_results:
                self.all_results['fcc_filtered_result'].to_excel(output_dir / (filename_prefix + '_fcc_result.xlsx'))
            print(f"Results saved to {output_dir}")
        except Exception as e:
            print(f"Failed to save results: {str(e)}")
            raise
    
    def get_model_info(self) -> Dict[str, Any]:
        """Get information about the loaded model"""
        if self.model is None:
            return {"status": "no_model_loaded"}
        
        info = {
            "model_type": self.parameters.model.get('model_type'),
            "model_path": self.parameters.model.get('model_path'),
            "solver": getattr(self.model, 'solver', 'unknown'),
            "tolerance": getattr(self.model, 'tolerance', 'unknown')
        }
        
        # Add model-specific information
        if hasattr(self.model, 'reactions'):
            info["num_reactions"] = len(self.model.reactions)
        if hasattr(self.model, 'metabolites'):
            info["num_metabolites"] = len(self.model.metabolites)
        if hasattr(self.model, 'genes'):
            info["num_genes"] = len(self.model.genes)
        
        return info
    
    def update_parameters(self, **kwargs) -> None:
        """Update workflow parameters"""
        model_updated = False

        for key, value in kwargs.items():
            if key in ('model', 'strain', 'algorithm') and isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    getattr(self.parameters, key)[sub_key] = sub_value
                    print(f"Updated parameter {key}.{sub_key} = {sub_value}")
                continue
            elif isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    if sub_key in self.parameters.model:
                        self.parameters.model[sub_key] = sub_value
                        print(f"Updated parameter model.{sub_key} = {sub_value}")
                        if sub_key == 'model_path':
                            model_updated = True
                    elif sub_key in self.parameters.strain:
                        self.parameters.strain[sub_key] = sub_value
                        print(f"Updated parameter strain.{sub_key} = {sub_value}")
                    elif sub_key in self.parameters.algorithm:
                        self.parameters.algorithm[sub_key] = sub_value
                        print(f"Updated parameter algorithm.{sub_key} = {sub_value}")
                    else:
                        print(f"Warning: Unknown parameter: {sub_key}")
                continue
            else:
                if key in self.parameters.model:
                    self.parameters.model[key] = value
                    print(f"Updated parameter model.{key} = {value}")
                    if key == 'model_path':
                        model_updated = True
                elif key in self.parameters.strain:
                    self.parameters.strain[key] = value
                    print(f"Updated parameter strain.{key} = {value}")
                elif key in self.parameters.algorithm:
                    self.parameters.algorithm[key] = value
                    print(f"Updated parameter algorithm.{key} = {value}")
                else:
                    print(f"Warning: Unknown parameter: {key}")
        
        if model_updated:
            self._validate_parameters()
            self._load_model()
    
    def get_results_summary(self) -> Dict[str, Any]:
        """Get a summary of the workflow results"""
        if not self.all_results:
            return {"status": "no_results"}
        
        summary = {
            "algorithm": self.parameters.algorithm.get('design_algorithm'),
            "target_reaction": self.parameters.strain.get('target_id'),
            "product": self.parameters.strain.get('product_name'),
            "status": "completed"
        }
        
        # Add algorithm-specific summary information
        if 'geneTable' in self.all_results:
            gene_table = self.all_results['geneTable']
            if isinstance(gene_table, pd.DataFrame):
                summary["num_targets"] = len(gene_table)
                summary["ko_targets"] = len(gene_table[gene_table.get('action', '') == 'KO'])
                summary["kd_targets"] = len(gene_table[gene_table.get('action', '') == 'KD'])
                summary["oe_targets"] = len(gene_table[gene_table.get('action', '') == 'OE'])
        
        return summary


