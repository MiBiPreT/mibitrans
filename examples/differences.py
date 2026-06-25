#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Sat Mar 14 12:09:46 2026.

@author: alraune
"""

import matplotlib.pyplot as plt
import numpy as np
import mibitrans as mbt


def mask(a, b, cutoff):
    """Set values in arrays a and b to NaN if both are below the provided cutoff."""
    if cutoff:
        mask = (a >= cutoff) | (b >= cutoff)
        return np.where(mask, a, np.nan), np.where(mask, b, np.nan)
    else:
        return a, b


def check_shape(a, b):
    """Check if a and b have the same shape."""
    if a.shape != b.shape:
        raise ValueError(f"cxyt_a and cxyt_b must have same shape, but have shapes a:{a.shape}, b:{b.shape}.")


def absolute_error(a, b):
    """Calculate the absolute error."""
    return abs(b - a)


def relative_error(a, b):
    """Calculate the relative error."""
    return abs(b - a) / a


def mean_absolute_difference(
    cxyt_a: np.ndarray, cxyt_b: np.ndarray, mean_axis: int | tuple = (0, 1, 2), concentration_cutoff: float = 1e-5
):
    """Calculate the mean absolute difference between two input arrays.

    Args:
        cxyt_a (np.ndarray) : Concentration array, must have same shape as cxyt_b.
        cxyt_b (np.ndarray) : Concentration array, must have same shape as cxyt_a.
        mean_axis (np.ndarray) : Over which axi to take the average difference. For example, `(1,2)` takes average over
            y and x dimensions, resulting in values for each time step. By default, averages over the entire array.
        concentration_cutoff (float) : If both cxyt have concentrations below this cutoff value, they will not be taken
            into account for calculating the average. If set to None, no cutoff value will be considered.
            Default is 1e-5.
    """
    check_shape(cxyt_a, cxyt_b)
    cxyt_a_masked, cxyt_b_masked = mask(cxyt_a, cxyt_b, concentration_cutoff)
    return np.nanmean(absolute_error(cxyt_a_masked, cxyt_b_masked), axis=mean_axis)


def mean_relative_difference(
    cxyt_a: np.ndarray, cxyt_b: np.ndarray, mean_axis: int | tuple = (0, 1, 2), concentration_cutoff: float = 1e-5
):
    """Calculate the mean relative difference between two input arrays.

    Args:
        cxyt_a (np.ndarray) : Concentration array, must have same shape as cxyt_b.
        cxyt_b (np.ndarray) : Concentration array, must have same shape as cxyt_a.
        mean_axis (np.ndarray) : Over which axi to take the average difference. For example, `(1,2)` takes average over
            y and x dimensions, resulting in values for each time step. By default, averages over the entire array.
        concentration_cutoff (float) : If both cxyt have concentrations below this cutoff value, they will not be taken
            into account for calculating the average. If set to None, no cutoff value will be considered.
            Default is 1e-5.
    """
    check_shape(cxyt_a, cxyt_b)
    cxyt_a_masked, cxyt_b_masked = mask(cxyt_a, cxyt_b, concentration_cutoff)
    return np.nanmean(relative_error(cxyt_a_masked, cxyt_b_masked), axis=mean_axis)

def mean_sum_absolute_difference(
    cxyt_a: np.ndarray, cxyt_b: np.ndarray, concentration_cutoff: float = 1e-5
):
    """Calculate the mean absolute difference of the sum of differences in the y-direction between two input arrays.

    Args:
        cxyt_a (np.ndarray) : Concentration array, must have same shape as cxyt_b.
        cxyt_b (np.ndarray) : Concentration array, must have same shape as cxyt_a.
        concentration_cutoff (float) : If both cxyt have concentrations below this cutoff value, they will not be taken
            into account for calculating the average. If set to None, no cutoff value will be considered.
            Default is 1e-5.
    """
    check_shape(cxyt_a, cxyt_b)
    y_length = len(cxyt_a[0,:,0])
    cxyt_a_center, cxyt_b_center = cxyt_a[:,y_length//2,:], cxyt_b[:,y_length//2,:]
    mask = (cxyt_a_center >= concentration_cutoff) | (cxyt_b_center >= concentration_cutoff)
    a, b = np.where(mask[:,None,:], cxyt_a, np.nan), np.where(mask[:,None,:], cxyt_b, np.nan)
    error = absolute_error(a, b)
    sum_error_xt = np.sum(error, axis=1)
    mean_sum_error_t = np.nanmean(sum_error_xt, axis=1)
    mean_sum_error = np.nanmean(mean_sum_error_t)
    return mean_sum_error

def rmse(
    cxyt_a: np.ndarray,
    cxyt_b: np.ndarray,
    axis: int | tuple = (0, 1, 2),
    concentration_cutoff: float = 1e-5,
):
    """Calculate the root-mean-square error between two input arrays."""
    check_shape(cxyt_a, cxyt_b)
    cxyt_a_masked, cxyt_b_masked = mask(cxyt_a, cxyt_b, concentration_cutoff)
    return np.sqrt(np.nanmean((absolute_error(cxyt_a_masked, cxyt_b_masked)) ** 2, axis=axis))



def sensitivity_models(
        hydrological_parameters,
        attenuation_parameters,
        source_parameters,
        model_parameters,
        x_dispersivity,
        y_dispersivity,
        z_dispersivity,
):
    """Repeatedly run models with different dispersivities for sensitivity analysis."""
    list_mibitrans = [[[0]*len(x_dispersivity)
                        for _ in range(len(y_dispersivity))]
                        for _ in range(len(z_dispersivity))]
    list_anatrans = [[[0]*len(x_dispersivity)
                        for _ in range(len(y_dispersivity))]
                        for _ in range(len(z_dispersivity))]
    list_bioscreen = [[[0]*len(x_dispersivity)
                        for _ in range(len(y_dispersivity))]
                        for _ in range(len(z_dispersivity))]

    mibitrans_object = mbt.Mibitrans(hydrological_parameters,
                                     attenuation_parameters,
                                     source_parameters,
                                     model_parameters)

    anatrans_object = mbt.Anatrans(hydrological_parameters,
                                   attenuation_parameters,
                                   source_parameters,
                                   model_parameters)

    bioscreen_object = mbt.Bioscreen(hydrological_parameters,
                                     attenuation_parameters,
                                     source_parameters,
                                     model_parameters)

    velocity = hydrological_parameters.velocity
    porosity = hydrological_parameters.porosity
    diffusion = hydrological_parameters.diffusion

    for i in range(len(x_dispersivity)):
        for j in range(len(y_dispersivity)):
            for k in range(len(z_dispersivity)):
                hydro_dispersivity = mbt.HydrologicalParameters(
                    velocity=velocity,
                    porosity=porosity,
                    alpha_x=x_dispersivity[i],
                    alpha_y=y_dispersivity[j],
                    alpha_z=z_dispersivity[k],
                    diffusion=diffusion,
                )
                print(i,j,k) # Print to track progress
                mibitrans_object.hydrological_parameters = hydro_dispersivity
                list_mibitrans[k][j][i] = mibitrans_object.run()
                anatrans_object.hydrological_parameters = hydro_dispersivity
                list_anatrans[k][j][i] = anatrans_object.run()
                bioscreen_object.hydrological_parameters = hydro_dispersivity
                list_bioscreen[k][j][i] = bioscreen_object.run()

    return list_mibitrans, list_anatrans, list_bioscreen

def comparison_plot(
    list_mibitrans, list_compare, x_dispersivity, y_dispersivity, z_dispersivity,
    time_index=None,
    difference_method="absolute",
    cutoff=1e-5,
    relative_concentration=True,
    relative_diff_diff=False,
    as_percentage=False,
    colorbar_label=None,
):
    """Plot difference between models, averaged over model domain, for all combinations of dispersivity."""
    ##### Plotting preferences #####
    value_decimals = 4
    figuresize = (7,5)
    colormap = "viridis"
    second_y_axis_location = -0.24
    x_label = r"Longitudinal dispersivity ($\alpha_L$) [m]"
    first_y_label = r"Transverse vertical dispersivity $\alpha_V$ [m]"
    second_y_label = r"Transverse horizontal dispersivity $\alpha_T$ [m]"
    #################################
    if time_index is None:
        time_index =  len(list_mibitrans[0][0][0].t)

    if relative_concentration:
        c_mode = "relative_cxyt"
    else:
        c_mode ="cxyt"

    difference_array = np.zeros((len(y_dispersivity) * len(z_dispersivity), len(x_dispersivity)))
    if difference_method not in ["relative", "rmse", "absolute", "sum_mean"]:
        print("Difference method not recognized, defaulting to absolute")

    for i in range(len(x_dispersivity)):
        for j in range(len(y_dispersivity)):
            for k in range(len(z_dispersivity)):
                if difference_method == "relative":
                    difference_array[len(z_dispersivity)*j+k,i] = mean_relative_difference(
                        getattr(list_mibitrans[k][j][i], c_mode)[:time_index,:,:],
                        getattr(list_compare[k][j][i], c_mode)[:time_index,:,:],
                        concentration_cutoff=cutoff
                    )
                elif difference_method == "rmse":
                    difference_array[len(z_dispersivity)*j+k,i] = rmse(
                        getattr(list_mibitrans[k][j][i], c_mode)[:time_index,:,:],
                        getattr(list_compare[k][j][i], c_mode)[:time_index,:,:],
                        concentration_cutoff=cutoff
                    )
                elif difference_method == "sum_mean":
                    difference_array[len(z_dispersivity)*j+k,i] = mean_sum_absolute_difference(
                        getattr(list_mibitrans[k][j][i], c_mode)[:time_index,:,:],
                        getattr(list_compare[k][j][i], c_mode)[:time_index,:,:],
                        concentration_cutoff=cutoff
                    )
                else:
                    difference_array[len(z_dispersivity)*j+k,i] = mean_absolute_difference(
                        getattr(list_mibitrans[k][j][i], c_mode)[:time_index,:,:],
                        getattr(list_compare[k][j][i], c_mode)[:time_index,:,:],
                        concentration_cutoff=cutoff
                    )

    if as_percentage:
        difference_array = difference_array * 100
    if relative_diff_diff:
        difference_array = difference_array / np.max(difference_array)

    z_ticks = []
    for i in range(len(y_dispersivity)):
        for j in range(len(z_dispersivity)):
            z_ticks.append(f"{z_dispersivity[j]}")

    fig, ax = plt.subplots(figsize=figuresize)
    imshow_aspect = str(1/len(z_dispersivity))
    im = ax.imshow(
        difference_array,
        cmap=colormap,
        aspect=imshow_aspect
    )
    ax2 = ax.secondary_yaxis(second_y_axis_location)
    middle_z_dispersivity_axis = len(z_dispersivity)/2-0.5

    second_y_tick_location = np.linspace(
        middle_z_dispersivity_axis,
        (len(y_dispersivity)-1)*len(z_dispersivity) + middle_z_dispersivity_axis,
        len(y_dispersivity))

    ax.set_xticks(range(len(x_dispersivity)), labels=x_dispersivity)
    ax.set_yticks(range(len(z_ticks)), labels=z_ticks)
    ax2.set_yticks(second_y_tick_location, labels=y_dispersivity)

    ax.hlines(second_y_tick_location[:-1]+len(z_dispersivity)/2,
              -0.5,
              len(x_dispersivity)-0.5,
              color="white",
              lw=2)
    ax.vlines(
        np.arange(0.5,len(x_dispersivity)-1,1),
        -0.5, len(y_dispersivity)*len(z_dispersivity)-0.5,
        color="white",
        lw=2
    )

    for i in range(len(x_dispersivity)):
        for j in range(len(y_dispersivity)*len(z_dispersivity)):
            ax.text(
                i, j, np.round(difference_array[j, i], decimals=value_decimals),
                ha="center",
                va="center",
                color="w"
            )

    ax.set_xlabel(x_label)
    ax.set_ylabel(first_y_label)
    ax2.set_ylabel(second_y_label)
    fig.colorbar(im, label=colorbar_label)
