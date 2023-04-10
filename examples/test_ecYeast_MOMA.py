import sys
import os
os.chdir('/Users/xluhon/Documents/GitHub/ETFLdesigner')
sys.path.append(r"/Users/xluhon/Documents/GitHub/ETFLdesigner/ETFLdesigner")

# load packages
from cobra.io import load_matlab_model, read_sbml_model
from pytfa.optim.utils import symbol_sum
from cobra.util.solver import set_objective
import pandas as pd
from ETFLdesigner.simulation.ecYeastFlux import *

# heme model
ecYeast=read_sbml_model('examples/models/yeast/heme_ecYeastGEM.xml')
product_name='heme_a'
product_id='EX_heme_a'

model2 = ecYeast.copy()
model2.reactions.get_by_id("r_1714_REV").bounds = (0, 100)  # D-glucose exchange (reversible)
model2.objective = {model2.reactions.r_2111: 1}
fba_solution1 = model2.optimize()
model2.objective = {model2.reactions.EX_heme_a: 1}
#model2.objective = {model2.reactions.r_1899: 1}
fba_solution2 = model2.optimize()
result2 = productionEnvolpe(model=model2, biomassRxn='r_2111', targetRxn='EX_heme_a',substrateRxn='r_1714_REV', modelType = 'EC')


# plot
import matplotlib.pyplot as plt # For python 3.6: sudo apt-get install python3.6-tk   TO INSTALL Tkinter
x1 = result2['biomass']
y1l= result2['p_lb']
y1 = result2['p_ub']
plt.figure(figsize=(4,3))
plt.plot(x1,y1, 'k',  color='red')
plt.plot(x1,y1l, 'k',  color='red')
#plt.plot((0.32, 0.32), (0, 0.35518), 'k-', color='red')
plt.xlabel('growth(/h)')
plt.ylabel('production(mmol/gDW.h)')
plt.show()












# pipeline to run moma based on Cobrapy
# prepare the geneID
geneID = []
for x in ecYeast.genes:
    print(x.id)
    geneID.append(x.id)

# prepare ecModel for calculation
model2 = ecYeast.copy()
#model2 = ecYeastMinimalMedia(model2)
qs=3
model2.reactions.get_by_id("r_1714_REV").bounds = (0, qs)  # D-glucose exchange (reversible)
model2.objective = {model2.reactions.r_2111: 1}
fba_solution = model2.optimize()
# loop
geneID = geneID[0:20]
fold_p = []
for i, x in enumerate(geneID):
    print(i)
    x = geneID[i]
    gene_remove = []
    gene_remove.append(x)
    mutModel = getModelWithRemoveGene(model0=model2, gene_remove0=gene_remove)
    try:
        fba_solution1 = moma(mutModel, solution=fba_solution, linear=True)
        fold0 = fba_solution1.fluxes['EX_heme_a'] / (qs * 6 * 0.13 / 3)
        print(fold0)
        fold_p.append(fold0)
    except:
        fold0 = None
        fold_p.append(fold0)
moma_result = pd.DataFrame({'gene':geneID, 'fold_p': fold_p})
# it seems that the remove of target gene obtained from moma could lead to the non cell growth in ecModel
# should further check the result

