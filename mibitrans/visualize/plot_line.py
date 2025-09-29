"""Author: Jorrit Bakker.

Module plotting a 3D matrix of contaminant plume concentrations as a line.
"""

import matplotlib.pyplot as plt
import numpy as np
from mibitrans.data.check_input import _check_model_type
from mibitrans.data.check_input import _time_check
from mibitrans.data.check_input import _x_check
from mibitrans.data.check_input import _y_check
from mibitrans.transport.domenico import Domenico


def centerline(model, y_position=0, time=None, **kwargs):
    """Plot center of contaminant plume as a line, at a specified time and, optionally, y position.

    Args:
        model : Model object from mibitrans.transport.
        y_position : y-position across the plume (transverse horizontal direction) for the plot.
            By default, the center of the plume at y=0 is plotted.
        time (float): Point of time for the plot. By default, last point in time is plotted.
        **kwargs : Arguments to be passed to plt.plot().

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

def transverse(model, x_position, time=None, **kwargs):
    """Plot concentration distribution as a line horizontal transverse to the plume extent.

    Args:
        model : Model object from mibitrans.transport.
        x_position : x-position along the plume (longitudinal direction) for the plot.
        time (float): Point of time for the plot. By default, last point in time is plotted.
        **kwargs : Arguments to be passed to plt.plot().
    """
    _check_model_type(model, Domenico)
    t_pos = _time_check(model, time)
    x_pos = _x_check(model, x_position)
    plot_array = model.cxyt[t_pos, :, x_pos]

    plt.plot(model.y, plot_array, **kwargs)

    plt.ylim((0, np.max(plot_array) + 1 / 10 * np.max(plot_array)))
    plt.xlabel("y-position [m]")
    plt.ylabel(r"Concentration [g/$m^{3}$]")
    plt.grid(True)

def breakthrough(model, x_position, y_position=0, **kwargs):
    """Plot contaminant breakthrough curve at given x and y position in model domain.

    Args:
        model : Model object from mibitrans.transport.
        x_position : x-position along the plume (longitudinal direction).
        y_position : y-position across the plume (transverse horizontal direction).
            By default, at the center of the plume (at y=0).
        **kwargs : Arguments to be passed to plt.plot().
    """
    _check_model_type(model, Domenico)
    x_pos = _x_check(model, x_position)
    y_pos = _y_check(model, y_position)
    plot_array = model.cxyt[:, y_pos, x_pos]

    plt.plot(model.t, plot_array, **kwargs)

    plt.ylim(0, np.max(plot_array) + 1 / 10 * np.max(plot_array))
    plt.xlabel("Time [days]")
    plt.ylabel(r"Concentration [g/$m^{3}$]")
    plt.grid(True)