'''Basic tutorial for ETFL models.
There are 4 types of ETFL models:
cEFL: Classical ME model(include expression and stoichiometry constraints) with constant biomass composition.
        examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json
vEFL: Classical ME model(include expression and stoichiometry constraints) with variable biomass composition.
        examples/models/yeast/yeast8_vEFL_2584_enz_64_bins__20231221_104514.json
cETFL: Thermokinetic-constrained ME model(include expression, stoichiometry and thermokinetic constraints) with constant biomass composition.
        examples/models/yeast/SlackModel yeast8_cETFL_2584_enz_64_bins__20231221_084634.json
vETFL: Thermokinetic-constrained ME model(include expression, stoichiometry and thermokinetic constraints) with variable biomass composition
        examples/models/yeast/SlackModel yeast8_vETFL_2584_enz_64_bins__20231221_105544.json
'''
# 1.load and save model
from etfl.io.json import save_json_model, load_json_model
model=load_json_model('examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json',solver='optlang-gurobi')
save_json_model(model,'your_filepath')

# basic configuration for solver
from etfl.optim.config import standard_solver_config
standard_solver_config(model)


# 2.basic simulation
from etfl.optim.utils import safe_optim

targetID='r_2111'
carbon_source='r_1714'
# set objective
model.objective = targetID
# FBA
solution = safe_optim(model)

# ppFBA-- Minimize total sum of enzyme usage FBA
from strainOptimizer.simulation.pprotFBA import ppFBA
solution = ppFBA(model,targetID,carbon_source,c_uptake=1,model_type='etfl')

# TODO: MOMA simulation

# customized objective setting
from pytfa.optim.utils import symbol_sum

# minimize the sum of specific reactionlist fluxes
object_rxnlist=[]
expr=symbol_sum([model.reactions.get_by_id(x).forward_variable+model.reactions.get_by_id(x).reverse_variable for x in object_rxnlist])
model.objective=expr
model.objective_direction='min'
model.optimize()


# 3. gene modification



# extract constraints and variables


# 调用蛋白质，核苷酸等比例
prot_ratio = model.interpolation_variable.prot_ggdw.variable.primal
mrna_ratio = model.interpolation_variable.mrna_ggdw.variable.primal
dna_ratio = model.interpolation_variable.dna_ggdw.variable.primal
lipid_ratio = model.interpolation_variable.lipid_ggdw.variable.primal
carbohydrate_ratio = model.interpolation_variable.carbohydrate_ggdw.variable.primal
ion_ratio = model.interpolation_variable.ion_ggdw.variable.primal




# 单敲见gene_essentiality_yETFL.py

# logger是一个记录器
model.logger.warning('一个warning')


# 重写了model这个class，涉及了复杂的多重继承







# 拿到变量，酶，mRNA，等
Enz_vars = model.get_variables_of_type(EnzymeVariable)




