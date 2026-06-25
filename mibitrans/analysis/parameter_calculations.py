"""Author: Jorrit Bakker.

Module containing various methods that takes a dictionary of parameters as input and calculates the proper values that
can be used in transport equations.
"""

import numpy as np


def calculate_biodegradation_capacity(electron_acceptors, utilization_factor):
    """Determine biodegradation capacity based on electron acceptor concentrations and utilization factor."""
    return np.sum(electron_acceptors.array / utilization_factor.array)


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


def _calculate_discharge_and_average_source_zone_concentration(
    hydrological_parameters, source_parameters, biodegradation_capacity=0
):
    """Calculate discharge through source zone and average source zone concentration, returned in respective order."""
    y_src = np.zeros(len(source_parameters.source_zone_boundary) + 1)
    y_src[1:] = source_parameters.source_zone_boundary
    c_src = source_parameters.source_zone_concentration
    Q = (
        hydrological_parameters.velocity
        * hydrological_parameters.porosity
        * source_parameters.depth
        * np.max(y_src)
        * 2
    )

    weighted_conc = np.zeros(len(source_parameters.source_zone_boundary))
    for i in range(len(source_parameters.source_zone_boundary)):
        weighted_conc[i] = (y_src[i + 1] - y_src[i]) * c_src[i]

    c0_avg = biodegradation_capacity + np.sum(weighted_conc) / np.max(y_src)

    return Q, c0_avg


def calculate_source_depletion(hydrological_parameters, source_parameters, biodegradation_capacity=0):
    """Calculate source depletion rate based on input parameters."""
    if source_parameters.total_mass != np.inf:
        Q, c0_avg = _calculate_discharge_and_average_source_zone_concentration(
            hydrological_parameters, source_parameters, biodegradation_capacity
        )
        k_source = Q * c0_avg / source_parameters.total_mass
    # If source mass is infinite, no source depletion is considered
    else:
        k_source = 0

    return k_source
