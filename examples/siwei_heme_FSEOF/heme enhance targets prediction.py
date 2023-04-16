import cobra
from cobra import Model, Reaction, Metabolite
import os
from os.path import join
import numpy as np
import pandas as pd
import copy
import re
model = cobra.io.load_matlab_model(r"C:\Users\Administrator\Desktop\Augest\y8.mat")


Ysx = 0.122*180/1000

def simulateGrowth(model, alpha):      # alpha是16等分的区间的某个点取值
    tmpmodel = model.copy()
    tmpmodel.reactions.get_by_id('r_1714').bounds = (0,0)
    tmpmodel.reactions.get_by_id('r_1714_REV').bounds = (0,1.0)  # 葡萄糖摄入1
    tmpmodel.reactions.get_by_id('r_1761').bounds = (0,1000)    # 文中封闭酒精吸收
    # model.reactions.get_by_id('r_1761_REV').bounds = (0,0)

    # max growth
    tmpmodel.objective = 'r_2111'
    sol = tmpmodel.optimize()

    # fix growth and max product
    tmpgrow = sol.objective_value * 0.999*alpha       # 令 tmgrow 等于最大biomass的0.999，这里我不理解为什么乘这么个alpha，有什么生理学意义？？？？？？
    tmpmodel.reactions.get_by_id('r_2111').bounds = (tmpgrow, tmpgrow)     # 固定biomass的值
    tmpmodel.objective = 'rEx'                                             # 换objective
    heme_max = tmpmodel.optimize()                                         # 这个optimize的结果名字叫heme_max
    # print(heme_max.objective_value)                                        # 看看结果

    # fix product and perform pFBA
    tmpmodel.reactions.get_by_id('rEx').bounds = (heme_max.objective_value * 0.999,
                                                  heme_max.objective_value * 0.999)  # 也固定heme_max的值
    sol_pfba = cobra.flux_analysis.pfba(tmpmodel)                                    # pFBA求出来
    return sol_pfba                                                                  # 输出pFBA结果

# 获得FSEOF的分区，如何取16个点（包括本身）
def div_(start_, end_, num):             #   这里的num要设置取的点数（不包括野生型本身的），其实就是几个间隔
    l = []
    interval = (end_ - start_)/num
    co = 0
    while co < num + 1:
        l.append(start_ + co * interval)
        co += 1
    return l          # 把r_2111的分布数列搞出来


def compare_substrate(model, Ysx):
    flux_WT = simulateGrowth(model, 1)        # 源代码，直接令Ysx是1，作为WT，俺也不知道为啥
    wtheme = Ysx / flux_WT.fluxes.loc['r_2111']     # pFBA里biomass那一行
    alpha = div_(wtheme / 2, wtheme * 2, 15)        # 这里最后会返回16个值
    v_matrix = []
    k_matrix = []
    f_gene = []
    for a in alpha:
        tmpflux = simulateGrowth(model, a)
        try:
            v_matrix.append(tmpflux.fluxes)
            k_matrix.append(np.asarray(tmpflux.fluxes)/     # 相同项数的数组除以数组，我没想到这个
                            np.asarray(flux_WT.fluxes))
        except AttributeError:                              # 有nan就会报错，不管它，承昱这么写就不报错可以出结果了
            v_matrix.append(tmpflux)
            k_matrix.append(tmpflux)
        # print(a)
        # print(alpha.index(a))


    return v_matrix, k_matrix, f_gene, alpha

v_matrix, k_matrix, f_gene, alpha = compare_substrate(model, Ysx)


all_flux = pd.DataFrame(v_matrix)   # 生成dataframe
all_flux = all_flux.T

all_k = pd.DataFrame(k_matrix)      # 生成dataframe
all_k = all_k.T
all_k.index = all_flux.index

# delete rxn without genes
f1 = list()
for h in all_flux.index:
    if bool(model.reactions.get_by_id(h).gene_reaction_rule == ''):
        f1.append(h)
    else:
        pass
all_flux.drop(index = f1, inplace = True)
all_k.drop(index = f1, inplace = True)

# delete all nan
all_flux.dropna(axis=0, how='all', inplace=True)
all_k.dropna(axis=0, how='all', inplace=True)

# nan --> 1
all_flux.fillna(1, inplace=True)
all_k.fillna(1, inplace=True)

# inf --> 1000
all_k.replace([np.inf, -np.inf], 1000, inplace=True)

# delete inconsistent
f5 = []                     # 这个过滤条件是承昱写的，hhhhh
for i in all_k.index:
    big = any(all_k.loc[i,:] > 1)       # 既有大于1的又有小于1的就不行
    small = any(all_k.loc[i,:] < 1)
    if big == small == True:
        f5.append(i)
all_k.drop(index = f5, inplace = True)

# order k(descend)
orderd_k = copy.deepcopy(all_k)
orderd_k['meank'] = orderd_k.mean(axis =1)
orderd_k = orderd_k.sort_values(by=["meank"], ascending=False)


# get gene，这部分我写的乱七八糟，直接照抄承昱写的理解了，循环里套循环，还要加上append，我还写不好
cons_g = {}  # 用一个字典，6
for ge in model.genes:
    # 取基因的k值
    ge_rxn = [r.id for r in ge.reactions if r.id in orderd_k.index]       # 用来储存对应orderd_k的index中的rxns
    ge_k = [orderd_k.loc[r.id, 'meank'] for r in ge.reactions if r.id in orderd_k.index]     # 拿出这个反应对应的k值
    # 看一致性
    if ge_k != []:    # 选择k值不是空的rxns
        gk = np.asarray(ge_k)    # 以数组形式提取出这个反应对应的k值，可以是1个或者好几个
        one = np.array(len(gk))  # 看这个反应，一共有几个相互对应的[基因，平均k值们]
        if sum(gk >= 1) == one or sum(gk <= 1) == one:    # 这里开始就没什么道理了，看不明白
            cons_g[ge.id] = np.mean(ge_k)

# 删mean到1
alpha_m = np.mean(alpha)
cons_g_f = {}
for k, v in cons_g.items():
    #if (v > 1e-3):   #文章源代码里没有这一步的
        if (v < alpha_m - 1e-3) or (v > 1 + 1e-3):
            cons_g_f[k] = v

# 排序
cons_g_f = pd.Series(cons_g_f)
cons_g_f.sort_values(ascending=False, inplace=True)

cons_g_f.to_excel(r"C:\Users\Administrator\Desktop\heme production enhance targets.xlsx")

