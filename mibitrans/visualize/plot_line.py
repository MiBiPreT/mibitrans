"""Author: Jorrit Bakker.

Module plotting a 3D matrix of contaminant plume concentrations as a line.
"""

import matplotlib.pyplot as plt
import numpy as np
from mibitrans.data.check_input import _check_model_type
from mibitrans.data.check_input import _time_check
from mibitrans.data.check_input import _y_check
from mibitrans.transport.domenico import Domenico


def centerline(model, time=None, y_position=0, **kwargs):
    """Plot center of contaminant plume as a line, at a specified time and, optionally, y position.

    Args:
        model : Model object from mibitrans.transport.
        time (float): Point of time for the plot. By default, last point in time is plotted.
        y_position : y-position across the plume (transverse horizontal direction) for the plot.
            By default, the center of the plume at y=0 is plotted.
        **kwargs : Arguments to be passed to plt.plot().

    Returns a line plot of the input plume as object.
    """
    _check_model_type(model, Domenico)
    t_pos = _time_check(model, time)
    y_pos = _y_check(model, y_position)
    plot_array = model.cxyt[t_pos, y_pos, :]

    plt.plot(model.x, plot_array, **kwargs)

    plt.ylim((0, np.max(plot_array) + 1 / 10 * np.max(plot_array)))
    plt.xlabel("Distance from source [m]")
    plt.ylabel(r"Concentration [g/$m^{3}$]")
    plt.grid(True)
