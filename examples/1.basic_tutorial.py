"""
1.  Some information
firstly: conda activate etfl
a)  Please run etfl models in the existing env named “etfl”.
b)  From now, only cplex is supported, so please double check your solver is cplex, not gurobi.
i.  for example: “/home/yeast/etfl/yetfl/code/benchmark_yeast_models” line 27, 28, 202
c)  All yeast model files with full functions are stored in “/home/yeast/etfl/yetfl/models”.
d)  The solving time of yeast etfl models various when single FBA is performed.
i.  cEFL/vEFL: 10min
ii. cETFL: 45-60 min
iii.    vETFL: 2-8 h
e)  For testing or cognition, please use the ecoli core model in “/home/yeast/etfl/etfl/tutorials/test_small.py”
"""


#2.  Usual codes

# load and save model
from etfl.io.json import save_json_model, load_json_model
load_json_model(path, solver=solver)
save_json_model(path)

# 调用蛋白质，核苷酸等比例
prot_ratio = model.interpolation_variable.prot_ggdw.variable.primal
mrna_ratio = model.interpolation_variable.mrna_ggdw.variable.primal
dna_ratio = model.interpolation_variable.dna_ggdw.variable.primal
lipid_ratio = model.interpolation_variable.lipid_ggdw.variable.primal
carbohydrate_ratio = model.interpolation_variable.carbohydrate_ggdw.variable.primal
ion_ratio = model.interpolation_variable.ion_ggdw.variable.primal

# 改solver
solver = 'optlang-gurobi'
solver = 'optlang-cplex'
model.solver = solver

# 最优化
standard_solver_config(model)
model.optimize()
# 或者用 
from etfl.optim.utils import safe_optim
out = safe_optim(model)

# 单敲见gene_essentiality_yETFL.py

# logger是一个记录器
model.logger.warning('一个warning')

# etfl model 有一个内置的solution,而GEM是没有的
model.solution

# 重写了model这个class，涉及了复杂的多重继承

# minimize substrate uptake
model.objective = symbol_sum([model.reactions.get_by_id(x).reverse_variable \
                              for x in anyUptake])
model.objective_direction = 'min'


# minimize total sum of fluxes
model.objective = symbol_sum([model.reactions.get_by_id(x.id).forward_variable + \
                              model.reactions.get_by_id(x.id).reverse_variable \
                              for x in cobra_model.reactions \
                              if x.id != 'r_4050']) 


# minimize enzyme usage i.e. max dummy enzyme
obj_expr = symbol_sum([model.enzymes.dummy_enzyme.variable])
set_objective(model,obj_expr)
model.objective_direction = 'max'

model.optimize()

# 拿到变量，酶，mRNA，等
Enz_vars = model.get_variables_of_type(EnzymeVariable)

# test config
from etfl.optim.config import standard_solver_config
standard_solver_config(model)
gene_ko_config(model)
growth_uptake_config(model)


