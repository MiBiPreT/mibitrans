"""Author: Jorrit Bakker.

Module including various methods to visualize (input) parameter conditions, intended to only be called internally.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mibitrans.analysis.parameter_calculations import calculate_biodegradation_capacity
from mibitrans.analysis.parameter_calculations import calculate_source_depletion
from mibitrans.data.check_input import check_instant_reaction_acceptor_input
from mibitrans.data.parameter_information import UtilizationFactor


def source_zone(source_parameters):
    """Visualize source zone conditions."""
    source_y = source_parameters.source_zone_boundary
    source_c = source_parameters.source_zone_concentration

    y_discretization = np.linspace(-source_y[-1] - source_y[-1] / 10, source_y[-1] + source_y[-1] / 10, 10000)
    c_values = np.zeros(len(y_discretization))
    for i, y in enumerate(source_y[::-1]):
        c_values = np.where((y_discretization <= y) & (y_discretization >= -y), source_c[-(i + 1)], c_values)

    indexer = np.linspace(1, 0.3, len(source_y))
    colormap = matplotlib.colormaps["YlGnBu"]

    for i, y in enumerate(source_y):
        if i == 0:
            plt.fill_betweenx(
                y=y_discretization,
                x1=c_values,
                where=(y_discretization <= y) & (y_discretization >= -y),
                color=colormap(indexer[i]),
                zorder=len(source_y) + 2,
            )
        else:
            plt.fill_betweenx(
                y=y_discretization,
                x1=c_values,
                # Boolean array for domain of source zone i
                where=((y_discretization <= y) & (y_discretization > source_y[i - 1]))
                | ((y_discretization >= -y) & (y_discretization < -source_y[i - 1])),
                color=colormap(indexer[i]),
                zorder=len(source_y) + 2 - i,
            )

    plt.xlabel(r"Source zone concentration $g/m^3$")
    plt.ylabel("Source zone y-coordinate")
    plt.title("Concentration distribution in the source zone")


def source_depletion(
    hydrological_parameters,
    source_parameters,
    electron_acceptors=None,
    utilization_factor=UtilizationFactor(
        util_oxygen=3.14, util_nitrate=4.9, util_ferrous_iron=21.8, util_sulfate=4.7, util_methane=0.78
    ),
    **kwargs,
):
    """Visualize source depletion."""
    if electron_acceptors is not None:
        ea, uf = check_instant_reaction_acceptor_input(electron_acceptors, utilization_factor)
        biodegradation_capacity = calculate_biodegradation_capacity(ea, uf)
    else:
        biodegradation_capacity = 0

    k_source = calculate_source_depletion(hydrological_parameters, source_parameters, biodegradation_capacity)
    source_half_time = np.log(2) / k_source

    t = np.linspace(0, source_half_time * 5, 500)
    source_conc_time = np.exp(-k_source * t)
    plt.plot(t, source_conc_time, **kwargs)
    plt.xlabel("Time [days]")
    plt.ylabel(r"Relative source concentration [$C/C_0$]")


def model_grid(model_parameters):
    """Visualize the model grid."""
    return None
