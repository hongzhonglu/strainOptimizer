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
from ETFLdesigner.manipulation.mainFunction import *



# heme model
ecYeast = read_sbml_model('examples/models/yeast/heme_ecYeastGEM.xml')
easy_model = produceRxnList(ecYeast)


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

# test MOMA based on enzyme abundance adjustment
rxnID =[rxn.id for rxn in model2.reactions]
rxnID = [x for x in rxnID if "draw_prot_" in x]

# prepare the geneID
geneID = []
for x in ecYeast.genes:
    print(x.id)
    geneID.append(x.id)

# prepare ecModel for calculation
def RefSimulation(model0, D0, prxn, qp):
    """
    This funcion is used to simulate the chemostat growth of yeast
    Actually this function is general to solve the ecGEMs
    :param model0: a ecGEMs
    :param D0: a growth rate
    :return: solution of fluxes
    """
    growth = D0
    with model0:
        # model0 = ecYeastMinimalMedia(model0)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth, growth)
        model0.reactions.get_by_id(prxn).bounds = (qp, qp)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 1000)  # open the glucose
        model0.objective = {model0.reactions.r_1714_REV: -1}
        solution2 = model0.optimize()
        GR = solution2.fluxes["r_1714_REV"]  # get the glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (GR, GR * 1.001)
        model0.objective = {model0.reactions.prot_pool_exchange: -1}
        solution3 = model0.optimize()
    return solution3
# run MOMA of protein resource adjust
def symbol_sum_MOMA(variables):
    from sympy import Add
    return Add(*variables)

model2 = ecYeast.copy()
# set tolerance
model2.tolerance = 1e-9
model2.reactions.get_by_id("r_1714_REV").bounds = (0, 100)
fba_solution = RefSimulation(model0=model2, D0=0.3, prxn='EX_heme_a', qp=0.002)
qs = fba_solution['r_1714_REV']

# get the reference condition
uniprotID = 'P08417'
proID = 'draw_prot_' + uniprotID
target_pro = fba_solution[proID]




mutModel = model2.copy()
mutModel.reactions.get_by_id(proID).bounds = (target_pro*5, float('inf'))
test = mutModel.optimize()

# try to rewrite moma for ecModels specially!!
"""Provide minimization of metabolic adjustment (MOMA)."""
from typing import TYPE_CHECKING, Optional
from optlang.symbolics import Zero, add
from cobra.util import solver as sutil
from cobra.flux_analysis.parsimonious import pfba
from cobra.core import Model, Solution

mutModel.solver = sutil.choose_solver(mutModel, qp=True)
prob = mutModel.problem
v = prob.Variable("moma_old_objective")
c = prob.Constraint(
    mutModel.solver.objective.expression - v,
    lb=0.0,
    ub=0.0,
    name="moma_old_objective_constraint",
)
to_add = [v, c]
mutModel.objective = prob.Objective(Zero, direction="min", sloppy=True)
obj_vars = []
for r in mutModel.reactions:
    if 'draw_prot_' in r.id:
        flux = fba_solution.fluxes[r.id]
        if flux > 0:
            print(r.id)
            dist = prob.Variable("moma_dist_" + r.id)
            """
            const = prob.Constraint(
                r.flux_expression - dist,
                lb=flux,
                ub=flux,
                name="moma_constraint_" + r.id,
            )
            """
            const = prob.Constraint(
                (r.flux_expression) / (1 - dist),
                lb=flux,
                ub=flux,
                name="moma_constraint_" + r.id,
            )
            to_add.extend([dist, const])
            obj_vars.append(dist ** 2)

mutModel.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)
### mutModel.add_cons_vars(to_add)#????????????????




# solve
test2 = mutModel.optimize()
test2.fluxes['r_1714_REV']
test2.fluxes['r_2111']




"""
# generate objective function use a new way
ss2 = []
for id in rxnID:
    print(fba_solution[id])
    if fba_solution[id] > 0:
        #xx = (ecoli.enzymes.get_by_id(id).variable - out_new[id])**2
        xx = (1-(model2.reactions.get_by_id(id).flux_expression)/(fba_solution[id]))**2
        ss2.append(xx)
add_objective1 = symbol_sum_MOMA(ss2)
new_objective = model2.problem.Objective(add_objective1, direction='min')
"""



# generate objective function use an old way
add_objective0 = []
for id in rxnID:
    print(fba_solution[id])
    if fba_solution[id] > 0:
        newone ='(1 - ' + 'model2.reactions.' + id +'.flux_expression' + '/' + str(fba_solution[id]) + ')**2'
        add_objective0.append(newone)
add_objective1 = " + ".join(add_objective0)
print(add_objective1)
new_objective2 = model2.problem.Objective(((1 - model2.reactions.draw_prot_P00127.flux_expression/3.856234203865542e-05)**2 + (1 - model2.reactions.draw_prot_P00128.flux_expression/3.856234203865542e-05)**2 + (1 - model2.reactions.draw_prot_P00163.flux_expression/3.856234203865542e-05)**2 + (1 - model2.reactions.draw_prot_P00360.flux_expression/4.593393564814084e-05)**2 + (1 - model2.reactions.draw_prot_P00401.flux_expression/1.2232812323417449e-05)**2 + (1 - model2.reactions.draw_prot_P00410.flux_expression/1.2232812323417449e-05)**2 + (1 - model2.reactions.draw_prot_P00420.flux_expression/1.2232812323417449e-05)**2 + (1 - model2.reactions.draw_prot_P00425.flux_expression/1.2232812323417449e-05)**2 + (1 - model2.reactions.draw_prot_P00427.flux_expression/1.2232812323417449e-05)**2 + (1 - model2.reactions.draw_prot_P00498.flux_expression/2.5549095744e-06)**2 + (1 - model2.reactions.draw_prot_P00549.flux_expression/5.5001907612315925e-06)**2 + (1 - model2.reactions.draw_prot_P00560.flux_expression/2.482009478835841e-06)**2 + (1 - model2.reactions.draw_prot_P00815.flux_expression/2.0191891416192e-07)**2 + (1 - model2.reactions.draw_prot_P00830.flux_expression/1.423828758481343e-05)**2 + (1 - model2.reactions.draw_prot_P00854.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_P00856.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_P00890.flux_expression/4.744839629004036e-06)**2 + (1 - model2.reactions.draw_prot_P00899.flux_expression/5.966859558900253e-06)**2 + (1 - model2.reactions.draw_prot_P00912.flux_expression/5.015635523400212e-08)**2 + (1 - model2.reactions.draw_prot_P00924.flux_expression/5.791675150944627e-06)**2 + (1 - model2.reactions.draw_prot_P00927.flux_expression/8.81855393099999e-08)**2 + (1 - model2.reactions.draw_prot_P00931.flux_expression/4.461172203480189e-09)**2 + (1 - model2.reactions.draw_prot_P00937.flux_expression/6.278639277150267e-06)**2 + (1 - model2.reactions.draw_prot_P00942.flux_expression/6.953047565135013e-08)**2 + (1 - model2.reactions.draw_prot_P00958.flux_expression/1.2866321406e-05)**2 + (1 - model2.reactions.draw_prot_P03962.flux_expression/1.5682518191699975e-08)**2 + (1 - model2.reactions.draw_prot_P03965.flux_expression/9.646584446699995e-09)**2 + (1 - model2.reactions.draw_prot_P04037.flux_expression/1.2232812323417449e-05)**2 + (1 - model2.reactions.draw_prot_P04039.flux_expression/2.4465930054842172e-05)**2 + (1 - model2.reactions.draw_prot_P04046.flux_expression/1.531603543589949e-07)**2 + (1 - model2.reactions.draw_prot_P04076.flux_expression/5.706616325999999e-09)**2 + (1 - model2.reactions.draw_prot_P04161.flux_expression/1.135698864299962e-06)**2 + (1 - model2.reactions.draw_prot_P04173.flux_expression/3.4077061986e-06)**2 + (1 - model2.reactions.draw_prot_P04801.flux_expression/5.9983429896e-06)**2 + (1 - model2.reactions.draw_prot_P04802.flux_expression/4.6199463519e-06)**2 + (1 - model2.reactions.draw_prot_P05150.flux_expression/2.03904465e-07)**2 + (1 - model2.reactions.draw_prot_P05373.flux_expression/1.796029364e-06)**2 + (1 - model2.reactions.draw_prot_P05375.flux_expression/2.747946523649831e-07)**2 + (1 - model2.reactions.draw_prot_P05626.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_P05694.flux_expression/4.558257738099363e-06)**2 + (1 - model2.reactions.draw_prot_P06168.flux_expression/2.364887467866e-05)**2 + (1 - model2.reactions.draw_prot_P06169.flux_expression/7.082659755448331e-07)**2 + (1 - model2.reactions.draw_prot_P06174.flux_expression/6.243136330000001e-10)**2 + (1 - model2.reactions.draw_prot_P06197.flux_expression/1.439517408663068e-06)**2 + (1 - model2.reactions.draw_prot_P06208.flux_expression/2.2363141316999997e-06)**2 + (1 - model2.reactions.draw_prot_P06633.flux_expression/6.898355186399999e-09)**2 + (1 - model2.reactions.draw_prot_P06773.flux_expression/1.4285519999999998e-10)**2 + (1 - model2.reactions.draw_prot_P06785.flux_expression/6.24996e-09)**2 + (1 - model2.reactions.draw_prot_P07143.flux_expression/3.856234203865542e-05)**2 + (1 - model2.reactions.draw_prot_P07149.flux_expression/8.28502818731814e-06)**2 + (1 - model2.reactions.draw_prot_P07172.flux_expression/2.22526465128e-09)**2 + (1 - model2.reactions.draw_prot_P07244.flux_expression/2.433671087399919e-06)**2 + (1 - model2.reactions.draw_prot_P07245.flux_expression/1.4015307489904197e-06)**2 + (1 - model2.reactions.draw_prot_P07251.flux_expression/1.423828758481343e-05)**2 + (1 - model2.reactions.draw_prot_P07255.flux_expression/1.2232812323417449e-05)**2 + (1 - model2.reactions.draw_prot_P07256.flux_expression/3.856234203865542e-05)**2 + (1 - model2.reactions.draw_prot_P07257.flux_expression/3.856234203865542e-05)**2 + (1 - model2.reactions.draw_prot_P07258.flux_expression/9.646584446699995e-09)**2 + (1 - model2.reactions.draw_prot_P07259.flux_expression/3.1082974250399965e-08)**2 + (1 - model2.reactions.draw_prot_P07262.flux_expression/1.8787326146748387e-05)**2 + (1 - model2.reactions.draw_prot_P07263.flux_expression/1.72456396272e-07)**2 + (1 - model2.reactions.draw_prot_P07264.flux_expression/6.653542466999999e-06)**2 + (1 - model2.reactions.draw_prot_P07277.flux_expression/8.461045755515579e-07)**2 + (1 - model2.reactions.draw_prot_P07283.flux_expression/6.567377271000001e-05)**2 + (1 - model2.reactions.draw_prot_P07285.flux_expression/1.1933594530500506e-06)**2 + (1 - model2.reactions.draw_prot_P07342.flux_expression/8.496356816400001e-06)**2 + (1 - model2.reactions.draw_prot_P07702.flux_expression/7.141006858199999e-07)**2 + (1 - model2.reactions.draw_prot_P07806.flux_expression/5.5061940492e-07)**2 + (1 - model2.reactions.draw_prot_P07807.flux_expression/5.0128500239994685e-09)**2 + (1 - model2.reactions.draw_prot_P08067.flux_expression/3.856234203865542e-05)**2 + (1 - model2.reactions.draw_prot_P08417.flux_expression/5.764367633852525e-07)**2 + (1 - model2.reactions.draw_prot_P08456.flux_expression/1.227140390013406e-08)**2 + (1 - model2.reactions.draw_prot_P08524.flux_expression/1.123025273121176e-05)**2 + (1 - model2.reactions.draw_prot_P08525.flux_expression/3.856234203865542e-05)**2 + (1 - model2.reactions.draw_prot_P08566.flux_expression/1.7734342356900066e-06)**2 + (1 - model2.reactions.draw_prot_P09436.flux_expression/2.4937260906e-07)**2 + (1 - model2.reactions.draw_prot_P09440.flux_expression/1.1188913462618206e-07)**2 + (1 - model2.reactions.draw_prot_P09457.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_P09624.flux_expression/1.1876172448026112e-05)**2 + (1 - model2.reactions.draw_prot_P09950.flux_expression/6.75541316e-06)**2 + (1 - model2.reactions.draw_prot_P10174.flux_expression/1.2232812323417449e-05)**2 + (1 - model2.reactions.draw_prot_P10356.flux_expression/1.9336988532e-08)**2 + (1 - model2.reactions.draw_prot_P10614.flux_expression/8.265970826616518e-07)**2 + (1 - model2.reactions.draw_prot_P10659.flux_expression/1.8316873211267483e-05)**2 + (1 - model2.reactions.draw_prot_P10869.flux_expression/1.9922472323999984e-07)**2 + (1 - model2.reactions.draw_prot_P11154.flux_expression/3.632686750655997e-06)**2 + (1 - model2.reactions.draw_prot_P11353.flux_expression/3.2684902e-05)**2 + (1 - model2.reactions.draw_prot_P11412.flux_expression/6.964188718514324e-07)**2 + (1 - model2.reactions.draw_prot_P11986.flux_expression/2.17526291995658e-06)**2 + (1 - model2.reactions.draw_prot_P12684.flux_expression/4.0041145255124777e-07)**2 + (1 - model2.reactions.draw_prot_P12695.flux_expression/1.1253742428441415e-05)**2 + (1 - model2.reactions.draw_prot_P12709.flux_expression/5.785644046984267e-07)**2 + (1 - model2.reactions.draw_prot_P13188.flux_expression/3.4270487964e-06)**2 + (1 - model2.reactions.draw_prot_P13298.flux_expression/8.276822027999987e-07)**2 + (1 - model2.reactions.draw_prot_P13663.flux_expression/2.3826541019999983e-07)**2 + (1 - model2.reactions.draw_prot_P14020.flux_expression/7.641914882700001e-07)**2 + (1 - model2.reactions.draw_prot_P15019.flux_expression/6.1447621421086744e-06)**2 + (1 - model2.reactions.draw_prot_P15180.flux_expression/6.5160226962e-07)**2 + (1 - model2.reactions.draw_prot_P15454.flux_expression/1.2401442803999144e-09)**2 + (1 - model2.reactions.draw_prot_P15496.flux_expression/2.0657207841346275e-10)**2 + (1 - model2.reactions.draw_prot_P15624.flux_expression/2.3104218036000002e-08)**2 + (1 - model2.reactions.draw_prot_P15625.flux_expression/2.3104218036000002e-08)**2 + (1 - model2.reactions.draw_prot_P15700.flux_expression/1.892353908599994e-08)**2 + (1 - model2.reactions.draw_prot_P16120.flux_expression/1.1084266431599991e-08)**2 + (1 - model2.reactions.draw_prot_P16387.flux_expression/1.1253742428441415e-05)**2 + (1 - model2.reactions.draw_prot_P16451.flux_expression/1.1253742428441415e-05)**2 + (1 - model2.reactions.draw_prot_P16550.flux_expression/6.210246712500001e-08)**2 + (1 - model2.reactions.draw_prot_P16603.flux_expression/8.265970826616518e-07)**2 + (1 - model2.reactions.draw_prot_P16622.flux_expression/1.8521377790000002e-06)**2 + (1 - model2.reactions.draw_prot_P16862.flux_expression/1.7019490769492524e-06)**2 + (1 - model2.reactions.draw_prot_P17423.flux_expression/1.108426643159999e-06)**2 + (1 - model2.reactions.draw_prot_P17505.flux_expression/3.5211590288183464e-06)**2 + (1 - model2.reactions.draw_prot_P18408.flux_expression/1.70339115375e-06)**2 + (1 - model2.reactions.draw_prot_P18544.flux_expression/1.0787223371999998e-05)**2 + (1 - model2.reactions.draw_prot_P19097.flux_expression/8.28502818731814e-06)**2 + (1 - model2.reactions.draw_prot_P19262.flux_expression/3.197919616e-07)**2 + (1 - model2.reactions.draw_prot_P19414.flux_expression/4.678466420907178e-06)**2 + (1 - model2.reactions.draw_prot_P20049.flux_expression/4.035307157999999e-07)**2 + (1 - model2.reactions.draw_prot_P20051.flux_expression/3.6045366848999946e-06)**2 + (1 - model2.reactions.draw_prot_P20967.flux_expression/3.197919616e-07)**2 + (1 - model2.reactions.draw_prot_P21147.flux_expression/2.092872178150437e-05)**2 + (1 - model2.reactions.draw_prot_P21264.flux_expression/2.129423871599929e-07)**2 + (1 - model2.reactions.draw_prot_P21306.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_P21375.flux_expression/8.645898264873414e-07)**2 + (1 - model2.reactions.draw_prot_P21592.flux_expression/1.683752525e-07)**2 + (1 - model2.reactions.draw_prot_P21954.flux_expression/1.8208802558246948e-06)**2 + (1 - model2.reactions.draw_prot_P22133.flux_expression/9.989571533850052e-09)**2 + (1 - model2.reactions.draw_prot_P22289.flux_expression/3.856234203865542e-05)**2 + (1 - model2.reactions.draw_prot_P22768.flux_expression/9.691707131999999e-07)**2 + (1 - model2.reactions.draw_prot_P22803.flux_expression/1.70339115375e-06)**2 + (1 - model2.reactions.draw_prot_P23254.flux_expression/3.847341259380847e-06)**2 + (1 - model2.reactions.draw_prot_P23542.flux_expression/2.3057759961629897e-07)**2 + (1 - model2.reactions.draw_prot_P24521.flux_expression/5.44985002506991e-06)**2 + (1 - model2.reactions.draw_prot_P25340.flux_expression/4.0691940008530056e-06)**2 + (1 - model2.reactions.draw_prot_P25628.flux_expression/1.822394532565951e-10)**2 + (1 - model2.reactions.draw_prot_P26637.flux_expression/8.223836220599999e-07)**2 + (1 - model2.reactions.draw_prot_P27472.flux_expression/6.0869262042e-08)**2 + (1 - model2.reactions.draw_prot_P27616.flux_expression/1.1001900680999633e-07)**2 + (1 - model2.reactions.draw_prot_P28239.flux_expression/3.024259237891865e-06)**2 + (1 - model2.reactions.draw_prot_P28240.flux_expression/4.242223308763842e-05)**2 + (1 - model2.reactions.draw_prot_P28272.flux_expression/1.0516584131999983e-06)**2 + (1 - model2.reactions.draw_prot_P28777.flux_expression/3.2195636460000126e-05)**2 + (1 - model2.reactions.draw_prot_P28789.flux_expression/8.890133320000001e-06)**2 + (1 - model2.reactions.draw_prot_P29509.flux_expression/2.0906594167500002e-07)**2 + (1 - model2.reactions.draw_prot_P29704.flux_expression/1.7032964768979921e-06)**2 + (1 - model2.reactions.draw_prot_P29952.flux_expression/3.8427056613e-07)**2 + (1 - model2.reactions.draw_prot_P30902.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_P31116.flux_expression/2.0028295703999986e-06)**2 + (1 - model2.reactions.draw_prot_P31688.flux_expression/1.3381826832e-07)**2 + (1 - model2.reactions.draw_prot_P32178.flux_expression/6.799075393499995e-08)**2 + (1 - model2.reactions.draw_prot_P32179.flux_expression/5.419988257500001e-07)**2 + (1 - model2.reactions.draw_prot_P32263.flux_expression/2.4692599763999998e-08)**2 + (1 - model2.reactions.draw_prot_P32264.flux_expression/1.4881816421999998e-06)**2 + (1 - model2.reactions.draw_prot_P32288.flux_expression/1.5011875296119935e-06)**2 + (1 - model2.reactions.draw_prot_P32340.flux_expression/5.098074100152159e-06)**2 + (1 - model2.reactions.draw_prot_P32347.flux_expression/2.49917482e-05)**2 + (1 - model2.reactions.draw_prot_P32377.flux_expression/3.7815217820436986e-06)**2 + (1 - model2.reactions.draw_prot_P32449.flux_expression/2.295986114250009e-07)**2 + (1 - model2.reactions.draw_prot_P32452.flux_expression/7.211228387999992e-09)**2 + (1 - model2.reactions.draw_prot_P32462.flux_expression/4.340446620666699e-05)**2 + (1 - model2.reactions.draw_prot_P32473.flux_expression/1.1253742428441415e-05)**2 + (1 - model2.reactions.draw_prot_P32476.flux_expression/3.697979580078504e-05)**2 + (1 - model2.reactions.draw_prot_P32582.flux_expression/3.1940178e-08)**2 + (1 - model2.reactions.draw_prot_P32622.flux_expression/1.6925193000007828e-08)**2 + (1 - model2.reactions.draw_prot_P32799.flux_expression/1.2232812323417449e-05)**2 + (1 - model2.reactions.draw_prot_P32895.flux_expression/4.291208303999998e-05)**2 + (1 - model2.reactions.draw_prot_P33312.flux_expression/0.00012500100000000002)**2 + (1 - model2.reactions.draw_prot_P33401.flux_expression/5.311003728e-07)**2 + (1 - model2.reactions.draw_prot_P33734.flux_expression/4.998811705199999e-08)**2 + (1 - model2.reactions.draw_prot_P36010.flux_expression/7.015689445795175e-07)**2 + (1 - model2.reactions.draw_prot_P36013.flux_expression/8.802631724999918e-08)**2 + (1 - model2.reactions.draw_prot_P36148.flux_expression/5.765751995904765e-05)**2 + (1 - model2.reactions.draw_prot_P36421.flux_expression/2.3071014773999995e-06)**2 + (1 - model2.reactions.draw_prot_P37254.flux_expression/3.107867999980976e-07)**2 + (1 - model2.reactions.draw_prot_P37291.flux_expression/2.8639234147904026e-05)**2 + (1 - model2.reactions.draw_prot_P37292.flux_expression/2.8924468171763553e-05)**2 + (1 - model2.reactions.draw_prot_P37299.flux_expression/3.856234203865542e-05)**2 + (1 - model2.reactions.draw_prot_P38009.flux_expression/3.5851565879999336e-06)**2 + (1 - model2.reactions.draw_prot_P38066.flux_expression/1.9290000000000003e-10)**2 + (1 - model2.reactions.draw_prot_P38077.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_P38139.flux_expression/3.6190457836269996e-08)**2 + (1 - model2.reactions.draw_prot_P38145.flux_expression/6.1728000000000004e-06)**2 + (1 - model2.reactions.draw_prot_P38169.flux_expression/2.491856400000707e-07)**2 + (1 - model2.reactions.draw_prot_P38221.flux_expression/5.202674097671341e-07)**2 + (1 - model2.reactions.draw_prot_P38426.flux_expression/1.3381826832e-07)**2 + (1 - model2.reactions.draw_prot_P38604.flux_expression/6.894531340762506e-07)**2 + (1 - model2.reactions.draw_prot_P38620.flux_expression/4.291208303999998e-05)**2 + (1 - model2.reactions.draw_prot_P38625.flux_expression/9.301907165399921e-06)**2 + (1 - model2.reactions.draw_prot_P38627.flux_expression/2.575929252599986e-07)**2 + (1 - model2.reactions.draw_prot_P38635.flux_expression/3.223437624e-09)**2 + (1 - model2.reactions.draw_prot_P38705.flux_expression/4.946213002499999e-06)**2 + (1 - model2.reactions.draw_prot_P38707.flux_expression/1.0476910311e-06)**2 + (1 - model2.reactions.draw_prot_P38708.flux_expression/3.0600785141999996e-06)**2 + (1 - model2.reactions.draw_prot_P38709.flux_expression/1.22202413952e-06)**2 + (1 - model2.reactions.draw_prot_P38720.flux_expression/1.2106935181811595e-05)**2 + (1 - model2.reactions.draw_prot_P38795.flux_expression/7.903686300002242e-08)**2 + (1 - model2.reactions.draw_prot_P38840.flux_expression/5.059952763136884e-07)**2 + (1 - model2.reactions.draw_prot_P38858.flux_expression/4.406935194270472e-05)**2 + (1 - model2.reactions.draw_prot_P38891.flux_expression/2.561020488e-08)**2 + (1 - model2.reactions.draw_prot_P38913.flux_expression/1.00644e-06)**2 + (1 - model2.reactions.draw_prot_P38972.flux_expression/1.6808662096199438e-06)**2 + (1 - model2.reactions.draw_prot_P38998.flux_expression/1.5839006849999998e-06)**2 + (1 - model2.reactions.draw_prot_P38999.flux_expression/9.0851042472e-08)**2 + (1 - model2.reactions.draw_prot_P39006.flux_expression/6.316110161134323e-05)**2 + (1 - model2.reactions.draw_prot_P39522.flux_expression/3.136768632e-06)**2 + (1 - model2.reactions.draw_prot_P39692.flux_expression/3.821657265e-08)**2 + (1 - model2.reactions.draw_prot_P39726.flux_expression/3.0263805798469704e-07)**2 + (1 - model2.reactions.draw_prot_P39954.flux_expression/1.3517479759280098e-05)**2 + (1 - model2.reactions.draw_prot_P40012.flux_expression/1.26859026e-07)**2 + (1 - model2.reactions.draw_prot_P40086.flux_expression/2.17052553e-08)**2 + (1 - model2.reactions.draw_prot_P40495.flux_expression/1.7516708759999997e-06)**2 + (1 - model2.reactions.draw_prot_P40545.flux_expression/2.1557297872799998e-07)**2 + (1 - model2.reactions.draw_prot_P40825.flux_expression/3.978030096e-07)**2 + (1 - model2.reactions.draw_prot_P40988.flux_expression/4.6827523680000004e-11)**2 + (1 - model2.reactions.draw_prot_P41338.flux_expression/2.0588277845966008e-08)**2 + (1 - model2.reactions.draw_prot_P41939.flux_expression/1.020306900996796e-06)**2 + (1 - model2.reactions.draw_prot_P41940.flux_expression/2.7658274340000004e-07)**2 + (1 - model2.reactions.draw_prot_P43567.flux_expression/4.274477430100195e-06)**2 + (1 - model2.reactions.draw_prot_P43619.flux_expression/8.430501600002392e-08)**2 + (1 - model2.reactions.draw_prot_P46655.flux_expression/2.2590872831999997e-05)**2 + (1 - model2.reactions.draw_prot_P46969.flux_expression/2.3389097814617435e-06)**2 + (1 - model2.reactions.draw_prot_P47011.flux_expression/6.0869262042e-08)**2 + (1 - model2.reactions.draw_prot_P47096.flux_expression/2.023313100000574e-08)**2 + (1 - model2.reactions.draw_prot_P47119.flux_expression/1.1135879999999999e-09)**2 + (1 - model2.reactions.draw_prot_P47125.flux_expression/4.215250800001196e-08)**2 + (1 - model2.reactions.draw_prot_P47143.flux_expression/4.1876431095454783e-07)**2 + (1 - model2.reactions.draw_prot_P47169.flux_expression/3.821657265e-08)**2 + (1 - model2.reactions.draw_prot_P47176.flux_expression/7.2290277726e-08)**2 + (1 - model2.reactions.draw_prot_P48015.flux_expression/3.0263805798469704e-07)**2 + (1 - model2.reactions.draw_prot_P48445.flux_expression/2.0051271549073957e-10)**2 + (1 - model2.reactions.draw_prot_P49089.flux_expression/6.22451058e-06)**2 + (1 - model2.reactions.draw_prot_P49095.flux_expression/3.0263805798469704e-07)**2 + (1 - model2.reactions.draw_prot_P49367.flux_expression/2.3822295107999998e-05)**2 + (1 - model2.reactions.draw_prot_P50094.flux_expression/1.6259870231999862e-06)**2 + (1 - model2.reactions.draw_prot_P50113.flux_expression/7.141006858199999e-07)**2 + (1 - model2.reactions.draw_prot_P50861.flux_expression/7.2150000000000004e-09)**2 + (1 - model2.reactions.draw_prot_P51601.flux_expression/1.3546994999917077e-07)**2 + (1 - model2.reactions.draw_prot_P52290.flux_expression/7.893336000000056e-13)**2 + (1 - model2.reactions.draw_prot_P52867.flux_expression/7.64187587244e-07)**2 + (1 - model2.reactions.draw_prot_P52910.flux_expression/1.6458215579265017e-06)**2 + (1 - model2.reactions.draw_prot_P53045.flux_expression/1.216862927647611e-08)**2 + (1 - model2.reactions.draw_prot_P53090.flux_expression/5.264489648536883e-07)**2 + (1 - model2.reactions.draw_prot_P53128.flux_expression/9.9711223939922e-08)**2 + (1 - model2.reactions.draw_prot_P53199.flux_expression/4.715606902005443e-10)**2 + (1 - model2.reactions.draw_prot_P53204.flux_expression/2.0988846000005956e-07)**2 + (1 - model2.reactions.draw_prot_P53848.flux_expression/1.1290081799930891e-07)**2 + (1 - model2.reactions.draw_prot_P53852.flux_expression/6.667141334999999e-08)**2 + (1 - model2.reactions.draw_prot_P54115.flux_expression/4.6680166660032375e-06)**2 + (1 - model2.reactions.draw_prot_P54839.flux_expression/4.464929333015452e-05)**2 + (1 - model2.reactions.draw_prot_P54885.flux_expression/1.9291517693999998e-07)**2 + (1 - model2.reactions.draw_prot_P61829.flux_expression/0.00046274119799517485)**2 + (1 - model2.reactions.draw_prot_P80210.flux_expression/2.0756706515999545e-06)**2 + (1 - model2.reactions.draw_prot_P81449.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_P81450.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_P81451.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_Q00055.flux_expression/4.6840184924594636e-08)**2 + (1 - model2.reactions.draw_prot_Q00764.flux_expression/1.3381826832e-07)**2 + (1 - model2.reactions.draw_prot_Q00955.flux_expression/2.7617744986835327e-06)**2 + (1 - model2.reactions.draw_prot_Q01217.flux_expression/4.34293938e-06)**2 + (1 - model2.reactions.draw_prot_Q01519.flux_expression/1.2232812323417449e-05)**2 + (1 - model2.reactions.draw_prot_Q02196.flux_expression/1.3549648702500003e-06)**2 + (1 - model2.reactions.draw_prot_Q03266.flux_expression/1.914172799988283e-07)**2 + (1 - model2.reactions.draw_prot_Q04066.flux_expression/3.1417713000008913e-09)**2 + (1 - model2.reactions.draw_prot_Q04119.flux_expression/7.344025800021786e-12)**2 + (1 - model2.reactions.draw_prot_Q04396.flux_expression/6.746511003534311e-09)**2 + (1 - model2.reactions.draw_prot_Q04409.flux_expression/6.450319104651652e-07)**2 + (1 - model2.reactions.draw_prot_Q04728.flux_expression/9.530329162799999e-07)**2 + (1 - model2.reactions.draw_prot_Q04952.flux_expression/4.2782752965e-06)**2 + (1 - model2.reactions.draw_prot_Q05506.flux_expression/9.340780775999999e-07)**2 + (1 - model2.reactions.draw_prot_Q05533.flux_expression/1.359524642012758e-06)**2 + (1 - model2.reactions.draw_prot_Q05911.flux_expression/2.6578573421399286e-06)**2 + (1 - model2.reactions.draw_prot_Q05979.flux_expression/1.44523665000041e-07)**2 + (1 - model2.reactions.draw_prot_Q06405.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_Q06817.flux_expression/1.90032489438e-06)**2 + (1 - model2.reactions.draw_prot_Q07500.flux_expression/3.3855495259498505e-06)**2 + (1 - model2.reactions.draw_prot_Q08548.flux_expression/1.5682944513134657e-07)**2 + (1 - model2.reactions.draw_prot_Q08911.flux_expression/1.3720084592547078e-05)**2 + (1 - model2.reactions.draw_prot_Q12040.flux_expression/1.5991570239917072e-06)**2 + (1 - model2.reactions.draw_prot_Q12055.flux_expression/1.2009001389969379e-05)**2 + (1 - model2.reactions.draw_prot_Q12109.flux_expression/2.1106320092999997e-06)**2 + (1 - model2.reactions.draw_prot_Q12122.flux_expression/4.8263142293999997e-05)**2 + (1 - model2.reactions.draw_prot_Q12154.flux_expression/1.379463539991556e-12)**2 + (1 - model2.reactions.draw_prot_Q12165.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_Q12189.flux_expression/1.8235205709692184e-05)**2 + (1 - model2.reactions.draw_prot_Q12198.flux_expression/9.854407035e-08)**2 + (1 - model2.reactions.draw_prot_Q12233.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_Q12349.flux_expression/4.6274119799517485e-05)**2 + (1 - model2.reactions.draw_prot_Q12362.flux_expression/8.4174e-07)**2 + (1 - model2.reactions.draw_prot_Q12452.flux_expression/7.69989159951558e-06)**2 + (1 - model2.reactions.draw_prot_Q12676.flux_expression/4.351395599973364e-07)**2 + (1 - model2.reactions.draw_prot_Q99258.flux_expression/2.3809800000000002e-08)**2), direction='min', sloppy=True)

# Overexpression
# The following script can't be used if we adjust the abundance of proteins
# get the reference condition
uniprotID = 'P23254'
proID = 'draw_prot_' + uniprotID
target_pro = fba_solution[proID]
mutModel = model2.copy()
mutModel.reactions.get_by_id(proID).bounds = (target_pro*3, float('inf'))
test = mutModel.optimize()
# mutModel.objective = symbol_sum_MOMA(ss2)
mutModel.objective = new_objective2
#mutModel.reactions.get_by_id("r_1714_REV").bounds = (3.46, 3.46)
sol4 = mutModel.optimize()
fold0 = sol4.fluxes['EX_heme_a'] / 0.002
print(fold0)









# knock-out
# loop
geneID = geneID[0:5]
geneID = ['YBR263W', 'YCL040W', 'YLR446W', 'YMR241W', 'YMR246W', 'YBR221C', 'YAL044C', 'YDR148C', 'YER070W', 'YFR053C', 'YGL253W', 'YGR180C', 'YJL026W']
fold_p = []
for i, x in enumerate(geneID):
    print(i)
    x = geneID[i]
    x = 'YBR263W'
    gene_remove = []
    gene_remove.append(x)
    mutModel = getModelWithRemoveGene(model0=model2, gene_remove0=gene_remove)
    mutModel = model2.copy()
    try:
        mutModel.objective = new_objective2
        #mutModel.reactions.get_by_id("r_1714_REV").bounds = (3.46, 3.46)
        sol4 = mutModel.optimize()
        fold0 = sol4.fluxes['EX_heme_a'] / 0.002
        print(fold0)
        fold_p.append(fold0)
    except:
        fold0 = None
        fold_p.append(fold0)
moma_result = pd.DataFrame({'gene':geneID, 'fold_p': fold_p})
# it seems that the remove of target gene obtained from moma could lead to the non cell growth in ecModel
# should further check the result
