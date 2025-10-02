import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from mibitrans.data.check_input import _check_model_type
from mibitrans.data.check_input import _time_check
from mibitrans.transport.domenico import Domenico


def plume_2d(model, time=None, animate=False, **kwargs):
    """Plot contaminant plume as a 2D colormesh, at a specified time.

    Args:
        model : Model object from mibitrans.transport.
        time (float): Point of time for the plot. Will show the closest time step to given value.
            By default, last point in time is plotted.
        animate (bool): If True, animation of contaminant plume until given time is shown.
        **kwargs : Arguments to be passed to plt.pcolormesh().

    Returns a matrix plot of the input plume as object.
    """
    _check_model_type(model, Domenico)
    t_pos = _time_check(model, time)
    if not animate:
        plt.pcolormesh(model.x, model.y, model.cxyt[t_pos, :, :], **kwargs)
        plt.xlabel("Distance from source (m)")
        plt.ylabel("Distance from plume center (m)")
        plt.colorbar(label=r"Concentration (g/$m^{3}$)")
    else:
        fig, ax = plt.subplots()
        mesh = ax.pcolormesh(model.x, model.y, model.cxyt[0, :, :], vmin=0, vmax=np.max(model.cxyt), **kwargs)
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.set_label("Concentration (C/C0)")
        ax.set_xlabel("Distance from source (m)")
        ax.set_ylabel("Distance from plume center (m)")

        def update(frame):
            mesh.set_array(model.cxyt[frame, :, :])
            ax.set_title(f"Concentration distribution at t={model.t[frame]} days")
            return mesh

        ani = animation.FuncAnimation(fig=fig, func=update, frames=len(model.t))
        return ani


def plume_3d(model, time=None, animate=False, **kwargs):
    """Plot contaminant plume as a 3D surface, at a specified time.

    Args:
        model : Model object from mibitrans.transport.
        time (float): Point of time for the plot. Will show the closest time step to given value.
            By default, last point in time is plotted.
        animate (bool): If True, animation of contaminant plume until given time is shown.
        **kwargs : Arguments to be passed to plt.plot_surface().

    Returns:
        ax (matplotlib.axes._axes.Axes) : Returns matplotlib axes object of plume plot.
    """
    _check_model_type(model, Domenico)
    t_pos = _time_check(model, time)

    if not animate:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        ax.plot_surface(model.xxx[t_pos, :, :], model.yyy[t_pos, :, :], model.cxyt[t_pos, :, :], **kwargs)
        ax.view_init(elev=30, azim=310)
        ax.set_xlabel("Distance from source (m)")
        ax.set_ylabel("Distance from plume center (m)")
        ax.set_zlabel(r"Concentration [$g/m^{3}$]")
        return ax

    else:
        if "cmap" not in kwargs and "color" not in kwargs:
            kwargs["color"] = "tab:blue"
        model_max = np.max(model.cxyt)
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surface = ax.plot_surface(
            model.xxx[0, :, :],
            model.yyy[0, :, :],
            model.cxyt[0, :, :],
            vmin=0,
            vmax=model_max,
            **kwargs,
        )
        ax.set_xlabel("Distance from source (m)")
        ax.set_ylabel("Distance from plume center (m)")

        # plot_surface creates a static surface; need to create new plot every time step
        def update(frame):
            nonlocal surface
            surface.remove()
            surface = ax.plot_surface(
                model.xxx[frame, :, :],
                model.yyy[frame, :, :],
                model.cxyt[frame, :, :],
                vmin=0,
                vmax=model_max,
                **kwargs,
            )
            ax.set_title(f"Concentration distribution at t={model.t[frame]} days")
            return surface

        ani = animation.FuncAnimation(fig=fig, func=update, frames=len(model.t))
        return ani
