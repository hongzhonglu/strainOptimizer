# strainOptimizer — Examples & Tutorials

This folder contains self-contained examples and step-by-step tutorials showing how to use
**strainOptimizer** for phenotype simulation and computational strain design with
enzyme-constrained GEMs (ecGEMs) and ETFL models.

---

## Prerequisites

| Requirement | Details |
|---|---|
| **Python** | 3.9 or 3.10 (recommended via conda) |
| **LP solver** | Gurobi 10.0 (preferred) **or** CPLEX — a valid license is required for either. Free academic licenses are available from [Gurobi](https://www.gurobi.com/academia/academic-program-and-licenses/) and [IBM](https://www.ibm.com/academic/). |
| **Package** | Install strainOptimizer in editable mode from the repo root: `pip install -e .` |

Solver strings used in the examples: `"optlang-gurobi"` or `"optlang-cplex"`.

---

## Quick-start notebooks (recommended entry point)

Interactive Jupyter notebooks with inline outputs and explanatory text cells.
Run them from the `examples/` directory.

| Notebook | Model | Target | Algorithm | Description |
|---|---|---|---|---|
| [`1.ecFactory_design_example.ipynb`](1.ecFactory_design_example.ipynb) | ecGEM (yeast ecYeastGEM) | 2-phenylethanol | ecFactory | End-to-end demonstration of the ecFactory algorithm: loads an enzyme-constrained yeast model, runs enzyme usage variability analysis (EUVA), and returns a ranked gene target list with knockout / overexpression calls. The reference example for new users. |
| [`2.ecFSEOF_design_example.ipynb`](2.ecFSEOF_design_example.ipynb) | ecGEM (yeast ecYeastGEM) | 2-phenylethanol | ecFSEOF | Demonstrates the ecFSEOF variant of ecFactory. Flux-scanning is performed at enforced objective values across a growth range to score enzyme utilisation changes, producing k-scores and gene action predictions comparable to the ecFactory notebook. |
| [`2.ecFactory_sclareol_design_example.ipynb`](2.ecFactory_sclareol_design_example.ipynb) | ecGEM (yeast + sclareol pathway) | Sclareol | ecFactory | Applies ecFactory to a yeast ecGEM that has been extended with a heterologous sclareol biosynthesis pathway. Shows how to combine pathway addition (`0.add_sclareol_pathway.py`) with automated strain design. |

---

## Strain design scripts

Equivalent Python scripts for running strain design workflows from the command line or in an IDE.

| Script | Model type | Organism / target | Algorithm | Description |
|---|---|---|---|---|
| [`1.ecFactory_design_ecoli.py`](1.ecFactory_design_ecoli.py) | ecGEM (`eciML1515_batch`) | *E. coli* / L-Tryptophan | ecFactory | Predicts gene overexpression and knockout targets for tryptophan overproduction in *E. coli*. Uses the enzyme-constrained iML1515 model. Includes FCC (flux control coefficient) calculation. The canonical ecFactory script for a prokaryotic host. |
| [`1.ecFactory_etfl_test.py`](1.ecFactory_etfl_test.py) | ETFL (`yeast8_cEFL`) | *S. cerevisiae* / 2-phenylethanol | ecFactory | Runs ecFactory on a yeast ETFL model, which additionally captures transcription and translation constraints. Illustrates how to set solver tolerances and timeouts for the more complex ETFL formulation. |
| [`1.ecFactory_core_etfl_test.py`](1.ecFactory_core_etfl_test.py) | ETFL (yeast core) | *S. cerevisiae* | ecFactory | Lightweight test of ecFactory on a reduced-size ETFL model. Useful for verifying the installation and solver setup before running full-genome models. |
| [`iBridge.py`](iBridge.py) | GEM / ecGEM / ETFL | *S. cerevisiae* / spermidine | iBridge | Demonstrates the iBridge strain design algorithm as an alternative to ecFactory. Loops over multiple simulation methods (MOMA, MOPA, ppFBA) and model types, then evaluates predictions against published experimental gene-modification targets. |

---

## Model preparation utilities

Scripts that extend or modify a model before running a design workflow.

| Script | Description |
|---|---|
| [`0.add_sclareol_pathway.py`](0.add_sclareol_pathway.py) | Adds a three-reaction sclareol biosynthesis pathway (farnesyltranstransferase → labdenyl-PP synthase → sclareol synthase) to a yeast ecGEM. Demonstrates how to create new `Reaction` and `Metabolite` objects, attach them to an existing ecGEM, and verify the maximum theoretical production yield before saving the extended model. |
| [`0.build_heme_ETFL_model.py`](0.build_heme_ETFL_model.py) | Extends a yeast ecGEM with a heme *a* demand reaction to enable heme production simulations. Shows how to locate an existing metabolite by ID (with fallback for different naming conventions across model formats) and add a boundary reaction to expose a commodity metabolite as an objective. |

---

## Omics integration

| Script | Description |
|---|---|
| [`GIMME.py`](GIMME.py) | Integrates RNA-seq transcriptomics data into a yeast ecGEM using the GIMME algorithm to produce a context-specific, expression-consistent model. Requires the `troppo` package (`pip install troppo`). Reads sample-averaged FPKM values and applies an expression threshold to prune low-evidence reactions. |

---

## Visualisation

| Script | Description |
|---|---|
| [`plot_envelope.py`](plot_envelope.py) | Calculates and plots growth–production envelopes (phenotypic phase planes) for a set of experimentally reported gene modifications (knockouts and overexpressions). Compares wild-type and engineered strain envelopes for 2-phenylethanol and spermidine in a yeast ecGEM. Requires `matplotlib`. |

---

## ETFL model tutorials (`ETFL_tutorials/`)

A progressive series of six tutorials covering the ETFL model interface.
All tutorials use the pre-built yeast cEFL model (`yeast8_cEFL_2584_enz_64_bins__20231221_083715.json`).

> **What is cEFL / vEFL / cETFL / vETFL?**
> The naming encodes the constraint level:
> - **c/v** — constant or variable biomass composition
> - **EFL** — includes stoichiometry + expression (transcription/translation) constraints
> - **ETFL** — EFL + thermodynamic feasibility constraints

| Script | Description |
|---|---|
| [`1.basic_tutorial.py`](ETFL_tutorials/1.basic_tutorial.py) | **Load and run an ETFL model.** Covers loading a pre-built ETFL JSON model, configuring the solver, running FBA and ppFBA (protein-parsimonious FBA), extracting macromolecule allocation ratios (protein, mRNA, lipid, carbohydrate fractions in g/gDW), and accessing enzyme/mRNA decision variables. Start here for any ETFL work. |
| [`2.gene_modification.py`](ETFL_tutorials/2.gene_modification.py) | **Simulate gene knockouts, knockdowns, and overexpressions.** Shows the correct ETFL approach: constraining the *translation reaction* upper/lower bound rather than modifying gene objects directly. Includes a reusable `gene_knock_out_down_up()` helper with ratio-based scaling. |
| [`3.concentration_predict.py`](ETFL_tutorials/3.concentration_predict.py) | **Extract predicted enzyme and mRNA concentrations.** Demonstrates how to read enzyme abundance (mmol/gDW) and mRNA levels from the ETFL solution, and how to convert enzyme concentrations to individual peptide concentrations using `enzymes_to_peptides_conc()`. Also shows how to map genes to their catalytic reactions. |
| [`4.add_reaction.py`](ETFL_tutorials/4.add_reaction.py) | **Add a heterologous pathway to an ETFL model.** Adds a two-step itaconate biosynthesis pathway (cis-aconitate decarboxylase + transport) to the yeast cEFL model and computes maximum itaconate production under a minimum growth constraint. Illustrates working with new metabolites, reactions, and demand reactions inside the ETFL framework. |
| [`5.update_kcat.py`](ETFL_tutorials/5.update_kcat.py) | **Update enzyme turnover numbers (kcat) and rebuild enzymatic constraints.** Shows how to modify `kcat_fwd` / `kcat_bwd` on an existing enzyme object and then rebuild the `ForwardCatalyticConstraint` / `BackwardCatalyticConstraint` expressions so the model reflects the new kinetic parameters without full reconstruction. |
| [`6.change_total_protein_constraint.py`](ETFL_tutorials/6.change_total_protein_constraint.py) | **Modify the total proteome allocation constraint.** Uses `constrain_enzymes()` to change the upper bound on total enzyme usage (g/gDW), simulating proteome reallocation conditions such as nutrient limitation or heterologous protein burden. Works for both ETFL and ecGEM model types. |

---

## ecGEM model tutorials (`ecmodel_tutorials/`)

| Script | Description |
|---|---|
| [`add_new_enzyme_into_prot_pool.py`](ecmodel_tutorials/add_new_enzyme_into_prot_pool.py) | **Add a new enzyme to the ecGEM protein pool.** Demonstrates the standard ecGEM pattern for introducing a heterologous enzyme: creating an enzyme pseudo-metabolite, adding a `draw_protein` reaction that draws from the shared `prot_pool[c]` according to the enzyme's molecular weight, and verifying the connection. Required when extending an ecGEM with reactions catalysed by foreign enzymes. |

---

## Archived / advanced examples (`ARCHIVE/`)

The `ARCHIVE/` folder contains earlier development scripts retained for reference.
They may require manual path adjustments and are not guaranteed to work with the current API.

---

## Expected outputs

Running any strain design example produces Excel result files in `examples/results/`:

| File pattern | Contents |
|---|---|
| `*_level_1_result.xlsx` | All candidate genes passing the EUVA enzyme-usage filter |
| `*_level_2_result.xlsx` | Subset passing gene essentiality filtering |
| `*_level_3_result.xlsx` | Final minimal candidate set after flux control coefficient ranking |
| `*_fcc_result.xlsx` | Flux control coefficients for all candidate reactions |
| `*_combined_result.xlsx` | Full merged table across all levels with gene names and recommended actions |

Pre-computed result files for tryptophan (*E. coli*), 2-phenylethanol (yeast), and
3-hydroxypropionic acid (yeast, xylose) are included for reference.
