"""Author: Jorrit Bakker.

Module containing various methods that takes a dictionary of parameters as input and calculates the proper values that
can be used in transport equations.
"""

import numpy as np
from mibitrans.data.parameter_information import util_to_conc_name


def calculate_utilization(model):
    """Function that calculates relative use of electron acceptors in biodegradation of BTEX."""
    util_factor = model._utilization_factor.dictionary
    biodeg_array = np.zeros(len(list(util_factor.keys())))
    util_array = np.zeros(len(biodeg_array))

    for i, (key, value) in enumerate(util_factor.items()):
        biodeg_array[i] = getattr(model._electron_acceptors, util_to_conc_name[key]) / value
        util_array[i] = value

    biodegradation_capacity = np.sum(biodeg_array)
    fraction_total = biodeg_array / biodegradation_capacity
    mass_fraction = fraction_total * util_array

    return mass_fraction


def calculate_discharge_and_average_source_zone_concentration(model):
    """Calculate discharge through source zone and average source zone concentration, returned in respective order."""
    if model.mode == "instant_reaction":
        bc = model.biodegradation_capacity
    else:
        bc = 0
    y_src = np.zeros(len(model.source_parameters.source_zone_boundary) + 1)
    y_src[1:] = model.source_parameters.source_zone_boundary
    c_src = model.source_parameters.source_zone_concentration
    Q = (
        model.hydrological_parameters.velocity
        * model.hydrological_parameters.porosity
        * model.source_parameters.depth
        * np.max(y_src)
        * 2
    )

    weighted_conc = np.zeros(len(model.source_parameters.source_zone_boundary))
    for i in range(len(model.source_parameters.source_zone_boundary)):
        weighted_conc[i] = (y_src[i + 1] - y_src[i]) * c_src[i]

    c0_avg = bc + np.sum(weighted_conc) / np.max(y_src)

    return Q, c0_avg
