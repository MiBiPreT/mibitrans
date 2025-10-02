"""Author: Jorrit Bakker.

Module including various methods to visualize (input) parameter conditions, intended to only be called internally.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


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

    plt.figure(dpi=300)

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


def model_grid(model_parameters):
    """Visualize the model grid."""
    return None
