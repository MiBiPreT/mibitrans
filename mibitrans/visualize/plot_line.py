"""Author: Jorrit Bakker.

Module plotting a 3D matrix of contaminant plume concentrations as a line.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
import mibitrans
from mibitrans.data.check_input import check_model_type
from mibitrans.data.check_input import check_time_in_domain
from mibitrans.data.check_input import check_x_in_domain
from mibitrans.data.check_input import check_y_in_domain


def centerline(model, y_position=0, time=None, legend_names=None, animate=False, **kwargs):
    """Plot center of contaminant plume of one or multiple models as a line, at a specified time and y position.

    Args:
        model : Model object from mibitrans.transport, or list of model objects.
        y_position (float): y-position across the plume (transverse horizontal direction) for the plot.
            By default, the center of the plume at y=0 is plotted.
        time (float): Point of time for the plot. Will show the closest time step to given value.
            By default, last point in time is plotted.
        legend_names (str | list): List of legend names as strings, in the same order as given models.
            By default, no legend is shown.
        animate (bool): If True, animation of contaminant plume until given time is shown. If multiple models are given
            as input, dt should be the same for each one to ensure accurate animation. Default is False.
        **kwargs : Arguments to be passed to plt.plot().

    """
    if not isinstance(model, list):
        model = [model]
    if not isinstance(legend_names, list) and legend_names is not None:
        legend_names = [legend_names]

    plot_array_list = []
    # Checks for list model input: dt should be equal, time should be smaller than the smallest end time, y_position
    # should be inside narrowest domain boundaries
    for mod in model:
        check_model_type(mod, mibitrans.transport.model_parent.Transport3D)
        y_pos = check_y_in_domain(mod, y_position)
        t_pos = check_time_in_domain(mod, time)
        _run_model_if_model_has_not_ran(mod)

        if animate:
            plot_array_list.append(mod.cxyt[:, y_pos, :])
        else:
            plot_array_list.append(mod.cxyt[t_pos, y_pos, :])

    # Non-animated plot
    if not animate:
        for i, mod in enumerate(model):
            if legend_names is not None:
                plt.plot(mod.x, plot_array_list[i], label=legend_names[i], **kwargs)
            else:
                plt.plot(mod.x, plot_array_list[i], **kwargs)

        plt.ylim(bottom=0)
        plt.xlabel("Distance from source [m]")
        plt.ylabel(r"Concentration [g/$m^{3}$]")
        if legend_names is not None:
            plt.legend()

    # Animated plot
    else:
        fig, ax = plt.subplots()
        plot_bin = []
        for i, mod in enumerate(model):
            if legend_names is not None:
                line = ax.plot(mod.x, plot_array_list[i][0, :], label=legend_names[i])[0]
            else:
                line = ax.plot(mod.x, plot_array_list[i][0, :], **kwargs)[0]
            plot_bin.append(line)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("Distance from source [m]")
        ax.set_ylabel(r"Concentration [g/$m^{3}$]")
        if legend_names is not None:
            ax.legend()

        def update(frame):
            for i, mod in enumerate(model):
                plot_bin[i].set_xdata(mod.x)
                plot_bin[i].set_ydata(plot_array_list[i][frame, :])
                ax.set_title(f"Concentration distribution at t={mod.t[frame]} days")
            return plot_bin

        ani = animation.FuncAnimation(fig=fig, func=update, frames=t_pos + 1)
        return ani


def transverse(model, x_position, time=None, legend_names=None, animate=False, **kwargs):
    """Plot concentration distribution as a line horizontal transverse to the plume extent.

    Args:
        model : Model object from mibitrans.transport, or list of model objects.
        x_position : x-position along the plume (longitudinal direction) for the plot.
        time (float): Point of time for the plot. Will show the closest time step to given value.
            By default, last point in time is plotted.
        legend_names (str | list): List of legend names as strings, in the same order as given models.
            By default, no legend is shown.
        animate (bool): If True, animation of contaminant plume until given time is shown. If multiple models are given
            as input, dt should be the same for each one to ensure accurate animation. Default is False.
        **kwargs : Arguments to be passed to plt.plot().
    """
    if not isinstance(model, list):
        model = [model]
    if not isinstance(legend_names, list) and legend_names is not None:
        legend_names = [legend_names]

    plot_array_list = []
    # Checks for list model input: dt should be equal, time should be smaller than the smallest end time, y_position
    # should be inside narrowest domain boundaries
    for mod in model:
        check_model_type(mod, mibitrans.transport.model_parent.Transport3D)
        x_pos = check_x_in_domain(mod, x_position)
        t_pos = check_time_in_domain(mod, time)
        _run_model_if_model_has_not_ran(mod)

        if animate:
            plot_array_list.append(mod.cxyt[:, :, x_pos])
        else:
            plot_array_list.append(mod.cxyt[t_pos, :, x_pos])

    if not animate:
        for i, mod in enumerate(model):
            if legend_names is not None:
                plt.plot(mod.y, plot_array_list[i], label=legend_names[i], **kwargs)
            else:
                plt.plot(mod.y, plot_array_list[i], **kwargs)

        plt.ylim(bottom=0)
        plt.xlabel("y-position [m]")
        plt.ylabel(r"Concentration [g/$m^{3}$]")
        if legend_names is not None:
            plt.legend()
    else:
        fig, ax = plt.subplots()
        plot_bin = []
        max_conc = 0
        for i, mod in enumerate(model):
            if legend_names is not None:
                line = ax.plot(mod.y, plot_array_list[i][0, :], label=legend_names[i])[0]
            else:
                line = ax.plot(mod.y, plot_array_list[i][0, :], **kwargs)[0]
            if np.max(plot_array_list[i]) > max_conc:
                max_conc = np.max(plot_array_list[i])
            plot_bin.append(line)
        ax.set_ylim(bottom=0, top=max_conc + max_conc / 10)
        ax.set_xlabel("y-position [m]")
        ax.set_ylabel(r"Concentration [g/$m^{3}$]")
        if legend_names is not None:
            ax.legend()

        def update(frame):
            for i, mod in enumerate(model):
                plot_bin[i].set_xdata(mod.y)
                plot_bin[i].set_ydata(plot_array_list[i][frame, :])
                ax.set_title(f"Concentration distribution at t={mod.t[frame]} days")
            return plot_bin

        ani = animation.FuncAnimation(fig=fig, func=update, frames=t_pos + 1)
        return ani


def breakthrough(model, x_position, y_position=0, legend_names=None, animate=False, **kwargs):
    """Plot contaminant breakthrough curve at given x and y position in model domain.

    Args:
        model : Model object from mibitrans.transport, or list of model objects.
        x_position : x-position along the plume (longitudinal direction).
        y_position : y-position across the plume (transverse horizontal direction).
            By default, at the center of the plume (at y=0).
        legend_names (str | list): List of legend names as strings, in the same order as given models.
            By default, no legend is shown.
        animate (bool): If True, animation of contaminant plume until given time is shown. If multiple models are given
            as input, dt should be the same for each one to ensure accurate animation. Default is False.
        **kwargs : Arguments to be passed to plt.plot().
    """
    if not isinstance(model, list):
        model = [model]
    if not isinstance(legend_names, list) and legend_names is not None:
        legend_names = [legend_names]

    plot_array_list = []
    # Checks for list model input: dt should be equal, time should be smaller than the smallest end time, y_position
    # should be inside narrowest domain boundaries
    for mod in model:
        check_model_type(mod, mibitrans.transport.model_parent.Transport3D)
        x_pos = check_x_in_domain(mod, x_position)
        y_pos = check_y_in_domain(mod, y_position)
        _run_model_if_model_has_not_ran(mod)
        plot_array_list.append(mod.cxyt[:, y_pos, x_pos])

    # Non animated plot
    if not animate:
        for i, mod in enumerate(model):
            if legend_names is not None:
                plt.plot(mod.t, plot_array_list[i], label=legend_names[i], **kwargs)
            else:
                plt.plot(mod.t, plot_array_list[i], **kwargs)

        plt.ylim(bottom=0)
        plt.xlabel("Time [days]")
        plt.ylabel(r"Concentration [g/$m^{3}$]")
        if legend_names is not None:
            plt.legend()

    # Animated plot
    else:
        fig, ax = plt.subplots()
        plot_bin = []
        max_conc = 0
        max_time = 0
        for i, mod in enumerate(model):
            if legend_names is not None:
                line = ax.plot(mod.t[0], plot_array_list[i][0], label=legend_names[i])[0]
            else:
                line = ax.plot(mod.t[0], plot_array_list[i][0], **kwargs)[0]

            # As plot extent is decided by first initiation of plot, ensure that axis concentration and time limits
            # are corresponding with their maximum values
            if mod.t[-1] > max_time:
                max_time = mod.t[-1]
            if np.max(plot_array_list[i]) > max_conc:
                max_conc = np.max(plot_array_list[i])
            plot_bin.append(line)
        ax.set_xlim(right=max_time)
        ax.set_ylim(bottom=0, top=max_conc + max_conc / 10)
        ax.set_xlabel("Time [days]")
        ax.set_ylabel(r"Concentration [g/$m^{3}$]")
        if legend_names is not None:
            ax.legend()

        def update(frame):
            for i, mod in enumerate(model):
                plot_bin[i].set_xdata(mod.t[:frame])
                plot_bin[i].set_ydata(plot_array_list[i][:frame])
                ax.set_title(f"Breakthrough curve at t={mod.t[frame]} days")
            return plot_bin

        ani = animation.FuncAnimation(fig=fig, func=update, frames=len(model[0].t))
        return ani


def _run_model_if_model_has_not_ran(model):
    if not model.has_run:
        if model.verbose:
            print("Model has not run yet, calculating cxyt...")
        model.run()
