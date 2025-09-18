"""Author: Jorrit Bakker.

Module containing various methods that takes a dictionary of parameters as input and calculates the proper values that
can be used in transport equations.
"""

import numpy as np
from mibitrans.data.parameter_information import util_to_conc_name


def calculate_utilization(model):
    """Function that calculates relative use of electron acceptors in biodegradation of BTEX."""
    util_factor = model.deg_pars.utilization_factor.dictionary
    biodeg_array = np.zeros(len(list(util_factor.keys())))
    util_array = np.zeros(len(biodeg_array))

    for i, (key, value) in enumerate(util_factor.items()):
        biodeg_array[i] = getattr(model.deg_pars, util_to_conc_name[key]) / value
        util_array[i] = value

    biodegradation_capacity = np.sum(biodeg_array)
    fraction_total = biodeg_array / biodegradation_capacity
    mass_fraction = fraction_total * util_array

    return mass_fraction
