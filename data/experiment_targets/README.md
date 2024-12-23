# Experimental validation of target dataset
## 1. Data statistics
| Product | Targets number | productID | Type                   |Reference|
|:--------|----------------|---------|------------------------|---|
|2-phenylethanol| 11             | r_1589  | aromatic alcohol       |10.1101/2023.01.31.526512|
|spermidine| 22             | r_2051  | polyamine              |10.1038/s41929-021-00631-z|
|free fatty acids| 18             |     |         lipid               |10.1016/j.cell.2018.07.013|
| heme| 38             |     |         porphyrin derivatives               | 10.1073/pnas.2108245119 |
|ergothioneine| 12             |  Heterologous       | amino acid derivatives |10.1016/j.ymben.2022.01.012|
|sclareol| 23             |  Heterologous       | terpene                |10.1016/j.ymben.2022.11.002|

## 2. Usage
```python
from strainOptimizer.analysis import dataset

# load experimental dataset
exp_data=dataset.load_experiment_targets(product='2-phenylethanol')

# Calculate experimental consistency
consistency=dataset.calculate_exp_consistency(predict_result=predict_result,exp_data=exp_data)
```



by: Haoyu Wang