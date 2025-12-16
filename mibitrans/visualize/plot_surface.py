import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from mibitrans.data.check_input import check_model_type
from mibitrans.data.check_input import check_time_in_domain
from mibitrans.visualize.plot_line import _plot_title_generator
from mibitrans.visualize.plot_line import allowed_model_types

relative_conc_zlabel = r"Relative concentration ($C/C_0$)"
absolute_conc_zlabel = r"Concentration [g/$m^{3}$]"


def plume_2d(model, time=None, relative_concentration=False, animate=False, **kwargs):
    """Plot contaminant plume as a 2D colormesh, at a specified time.

    Args:
        model : Model object from mibitrans.transport.
        time (float): Point of time for the plot. Will show the closest time step to given value.
            By default, last point in time is plotted.
        relative_concentration (bool, optional) : If set to True, will plot concentrations relative to maximum source
            zone concentrations at t=0. By default, absolute concentrations are shown.
        animate (bool, optional): If True, animation of contaminant plume until given time is shown. Default is
            False.
        **kwargs : Arguments to be passed to plt.pcolormesh().

    Returns a matrix plot of the input plume as object.
    """
    check_model_type(model, allowed_model_types())
    t_pos = check_time_in_domain(model, time)
    if relative_concentration:
        model_concentration = model.relative_cxyt
        z_label = relative_conc_zlabel
    else:
        model_concentration = model.cxyt
        z_label = absolute_conc_zlabel
    # Non animated plot
    if not animate:
        plt.pcolormesh(model.x, model.y, model_concentration[t_pos, :, :], **kwargs)
        plt.xlabel("Distance from source (m)")
        plt.ylabel("Distance from plume center (m)")
        plt.colorbar(label=z_label)
        plot_title = _plot_title_generator("Plume", model, time=model.t[t_pos])
        plt.title(plot_title)

    # Animated plot
    else:
        fig, ax = plt.subplots()
        mesh = ax.pcolormesh(
            model.x, model.y, model_concentration[0, :, :], vmin=0, vmax=np.max(model_concentration), **kwargs
        )
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.set_label(z_label)
        ax.set_xlabel("Distance from source (m)")
        ax.set_ylabel("Distance from plume center (m)")

        def update(frame):
            mesh.set_array(model_concentration[frame, :, :])
            ax.set_title(f"Concentration distribution at t={model.t[frame]} days")
            return mesh

        ani = animation.FuncAnimation(fig=fig, func=update, frames=t_pos + 1)
        return ani


def plume_3d(model, time=None, relative_concentration=False, animate=False, **kwargs):
    """Plot contaminant plume as a 3D surface, at a specified time.

    Args:
        model : Model object from mibitrans.transport.
        time (float): Point of time for the plot. Will show the closest time step to given value.
            By default, last point in time is plotted.
        relative_concentration (bool, optional) : If set to True, will plot concentrations relative to maximum source
            zone concentrations at t=0. By default, absolute concentrations are shown.
        animate (bool, optional): If True, animation of contaminant plume until given time is shown. Default is
            False.
        **kwargs : Arguments to be passed to plt.plot_surface().

    Returns:
        ax (matplotlib.axes._axes.Axes) : Returns matplotlib axes object of plume plot.
    """
    check_model_type(model, allowed_model_types())
    t_pos = check_time_in_domain(model, time)
    if relative_concentration:
        model_concentration = model.relative_cxyt
        z_label = relative_conc_zlabel
    else:
        model_concentration = model.cxyt
        z_label = absolute_conc_zlabel

    # Non animated plot
    xxx = np.tile(model.x, (len(model.t), len(model.y), 1))
    yyy = np.tile(model.y[:, None], (len(model.t), 1, len(model.x)))
    if not animate:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        ax.plot_surface(xxx[t_pos, :, :], yyy[t_pos, :, :], model_concentration[t_pos, :, :], **kwargs)
        ax.view_init(elev=30, azim=310)
        ax.set_xlabel("Distance from source (m)")
        ax.set_ylabel("Distance from plume center (m)")
        ax.set_zlabel(z_label)
        plot_title = _plot_title_generator("Plume", model, time=model.t[t_pos])
        ax.set_title(plot_title)
        return ax

    # Animated plot
    else:
        if "cmap" not in kwargs and "color" not in kwargs:
            kwargs["color"] = "tab:blue"
        model_max = np.max(model_concentration)
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surface = ax.plot_surface(
            xxx[0, :, :],
            yyy[0, :, :],
            model_concentration[0, :, :],
            vmin=0,
            vmax=model_max,
            **kwargs,
        )
        ax.set_xlabel("Distance from source (m)")
        ax.set_ylabel("Distance from plume center (m)")
        ax.set_zlabel(z_label)
        ax.set_zlim(0, model_max)

        # plot_surface creates a static surface; need to create new plot every time step
        def update(frame):
            # nonlocal needed in order for the previous plot to be removed before new one is plotted
            nonlocal surface
            surface.remove()
            surface = ax.plot_surface(
                xxx[frame, :, :],
                yyy[frame, :, :],
                model_concentration[frame, :, :],
                vmin=0,
                vmax=model_max,
                **kwargs,
            )
            ax.set_title(f"Concentration distribution at t={model.t[frame]} days")
            return surface

        ani = animation.FuncAnimation(fig=fig, func=update, frames=t_pos + 1)
        return ani
