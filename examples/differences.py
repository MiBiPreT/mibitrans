#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Sat Mar 14 12:09:46 2026.

@author: alraune
"""

import numpy as np


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
