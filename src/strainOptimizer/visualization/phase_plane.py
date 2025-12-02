# -*- coding: utf-8 -*-
# date : 2024/5/9 
# author : wangh
from matplotlib import pyplot as plt
import numpy as np


def calculate_flux_space(model,x_axis,y_axis,points=20):
    ''' Calculate the flux space relationship between rxn x and y
    *param:
        model: cobra model
        x_axis: rxn id
        y_axis: rxn id
        points: number of points
    *return:
        flux_space_dict: dict,
            x_values: flux values for rxn x
            y_ubs: flux upper bounds for rxn y
            y_lbs: flux lower bounds for rxn y
            x_axis: rxn id for x
            y_axis: rxn id for y
    * example:
    flux_space=calculate_flux_space(model,x_axis='r_2111',y_axis='r_2154',points=20)

    '''
    model.objective = x_axis
    model.objective_direction = 'max'
    x_max=model.slim_optimize()
    model.objective_direction = 'min'
    x_min=model.slim_optimize()
    x_values = np.linspace(x_min,x_max,points)
    y_ubs = []
    y_lbs = []
    for x in x_values:
        model.reactions.get_by_id(x_axis).bounds = x, x
        model.objective = y_axis
        model.objective_direction = 'max'
        y_ub = model.slim_optimize()
        y_ubs.append(y_ub)
        model.objective_direction = 'min'
        y_lb = model.slim_optimize()
        y_lbs.append(y_lb)

    flux_space_dict={'x_values':x_values,'y_ubs':y_ubs,'y_lbs':y_lbs,'x_axis':x_axis,'y_axis':y_axis}
    return flux_space_dict

def plot_envelope(flux_spaceList,labels,show=True):
    '''Plot the phase plane for target 2 reactions.
    *param:
        flux_spaceList: list of flux_space_dict
        labels: list of labels
        show: bool, whether to show the plot

    *return:
        fig: figure object

    * example:
    flux_spaceList=[calculate_flux_space(model1,x_axis='r_2111',y_axis='r_2154',points=20),
            calculate_flux_space(model2,x_axis='r_2111',y_axis='r_2154',points=20)]
    labels=['wild-type','mutant']
    fig=plot_envelope(flux_spaceList,labels,show=True)
    '''
    # check if flux_spaceList and labels are the same length
    if len(flux_spaceList) != len(labels):
        raise ValueError('flux_spaceList and labels must be the same length')
    fig, ax = plt.subplots()
    for i in range(len(flux_spaceList)):
        flux_space = flux_spaceList[i]
        label=labels[i]
        # set the line width
        ax.fill_between(flux_space['x_values'], flux_space['y_ubs'], flux_space['y_lbs'], alpha=0.5,label=label,
                        linewidth=2,
                        # edgecolor='black'
                        )
    # set axis
    ax.set_xlabel(flux_space['x_axis'])
    ax.set_ylabel(flux_space['y_axis'])
    # set x,y limit as 0
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    # set legend
    ax.legend()
    if show:
        plt.show()
    return fig
