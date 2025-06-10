# newly added
import sys
import os

# sys.path.append('/Users/xluhon/Documents/GitHub/etfl')
# sys.path.append(r"/Users/xluhon/Documents/GitHub/strainOptimizer/")
from strainOptimizer.manipulation.model_process import *
from strainOptimizer.manipulation.protein_process import *

from collections import namedtuple
import pandas as pd
import numpy  as np

from etfl.io.json import load_json_model
from pytfa.optim.utils import symbol_sum
from cobra.util.solver import set_objective

solver = 'optlang-gurobi'

# ec_cobra.reactions.ATPM.lower_bound = 0
growth_reaction_id = 'r_4041'
yeast = load_json_model('examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json', solver=solver)

# initial test
yeast.reactions.r_1714.lower_bound = -10
yeast.reactions.r_1714.upper_bound = 0
yeast.objective = growth_reaction_id
yeast.objective.direction = 'max'
sol = yeast.optimize()

# input the experimental datasets
# This physiology datasets are from Carl, PNAS, 2021.
# r_2111 growth 0.424
# r_1714_REV glucose uptake 18.79
# r_1992_REV oxygen uptake 2.92
# r_1672 co2 production 24.77
# r_1761 ethanol production 26.97
# r_1634 acetate secretion 0.566
# Glycerol secretion rate 1.469
# Pyruvate 0.163
# Succinate 0.0365
# growth 0.424
# protein_ratio 0.483


yeast.reactions.get_by_id("r_1634").upper_bound = 0 # assume acetate is not produced!
yeast.reactions.get_by_id("r_2033").upper_bound = 0 # assume pyruvate is not produced!
yeast.reactions.get_by_id("r_1631").upper_bound = 0 # assume acetaldehyde is not produced!
yeast.reactions.get_by_id("r_1549").upper_bound = 0 # assume (R,R)-2,3-butanediol is not produced!
#yeast.reactions.get_by_id("prot_pool_exchange").upper_bound = 0.10366*1.125
#yeast.reactions.get_by_id("prot_pool_exchange").upper_bound = 0.10366*1.125


# next using the chemostat simulation
growth = 0.424
yeast.reactions.get_by_id("r_2111").bounds = (growth, growth)

# minimization glucose uptake rate
yeast.reactions.get_by_id("r_1714").bounds = (-1000, 1000)  # open the glucose
yeast.objective = {yeast.reactions.r_1714: 1}
solution2 = yeast.optimize()

GR = solution2.fluxes["r_1714"]  # get the glucose uptake rate
yeast.reactions.get_by_id("r_1714").bounds = (GR * 1.001, GR)

obj_expr = symbol_sum([yeast.enzymes.dummy_enzyme.variable])
set_objective(yeast, obj_expr)
yeast.objective_direction = 'max'
yeast.optimize()


# how to output the pipetide concentration
from etfl.analysis.utils import enzymes_to_peptides_conc
ret = dict()
for enz in yeast.enzymes:
    ret['EZ_' + enz.id] = yeast.solution.raw.loc['EZ_' + enz.id]
peptides_conc = enzymes_to_peptides_conc(yeast, ret)
pro_abundance = pd.Series(peptides_conc)
result = pd.DataFrame({'rxnID':pro_abundance.index, 'flux':pro_abundance.values})
result.to_excel("/Users/xluhon/Documents/GitHub/strainOptimizer/examples/data/predicted_protein_abundance.xlsx")




# input the proteomics under max growth rate
result=pd.read_excel("/Users/xluhon/Documents/GitHub/strainOptimizer/examples/data/predicted_protein_abundance.xlsx")
abundance_ex = pd.read_excel("/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/data/proteomics/data_PNAS_2021.xlsx")
abundance_ex['g/gDW'] =(abundance_ex['replicate 1 (g gDW-1)']+ abundance_ex['replicate 2 (g gDW-1)']+ abundance_ex['replicate 3 (g gDW-1)'])/3
abundance_ex=abundance_ex[['Symbol','g/gDW']]
abundance_ex.columns = ['gene','g/gDW']
abundance_ex1 = splitAbundance(pro_df=abundance_ex)
# change the unit from g/gDW as mmol/gDW
# input the molecular weight

mw = pd.read_csv("/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000
abundance_ex1["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abundance_ex1["gene"])
abundance_ex_check = abundance_ex1[abundance_ex1["MW_Kda"].isna()]
abundance_ex1=abundance_ex1[~abundance_ex1["MW_Kda"].isna()]
abundance_ex1["mmol/gDW"] = abundance_ex1["g/gDW"]/abundance_ex1["MW_Kda"]# #mmol/g biomass



result['pro_measured'] = singleMapping(abundance_ex1["mmol/gDW"],abundance_ex1["gene"],result['rxnID'])
result.to_excel("/Users/xluhon/Documents/GitHub/strainOptimizer/examples/result/predicted_and_measured_proteomics.xlsx")
result=result[~result["pro_measured"].isna()]



# plot
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

plt.figure()
sns.regplot(x=result['pro_measured'], y=result['flux'], fit_reg=False)
sns.regplot(x=np.log10(result['pro_measured']), y=np.log10(result['flux']), fit_reg=False)
plt.xlim(-10, -3)
plt.ylim(-10, -3)
plt.xlabel("log10(Measured_protein_level)")
plt.ylabel("log10(Predicted_protein_usage)")


# remove the proteins with zero
result1 = result[result['flux'] > 0]
result1 = result1[result1['pro_measured'] > 0]
corr, ss = pearsonr(np.log10(result1['pro_measured']), np.log10(result1['flux']))
print("Correlation coefficient:", corr)
print("Correlation p_value:", ss)

