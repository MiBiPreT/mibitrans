import matplotlib.pyplot as plt
from mibitrans.data.check_input import _check_model_type
from mibitrans.data.check_input import _time_check
from mibitrans.transport.domenico import Domenico


def plume_2d(model, time=None, **kwargs):
    """Plot contaminant plume as a 2D colormesh, at a specified time.

    Args:
        model : Model object from mibitrans.transport.
        time (float): Point of time for the plot. By default, last point in time is plotted.
        **kwargs : Arguments to be passed to plt.pcolormesh().

    Returns a matrix plot of the input plume as object.
    """
    _check_model_type(model, Domenico)
    t_pos = _time_check(model, time)
    plt.pcolormesh(model.x, model.y, model.cxyt[t_pos, :, :], **kwargs)

    plt.xlabel("Distance from source (m)")
    plt.ylabel("Distance from plume center (m)")
    plt.colorbar(label=r"Concentration (g/$m^{3}$)")


def plume_3d(model, time=None, **kwargs):
    """Plot contaminant plume as a 3D surface, at a specified time.

    Args:
        model : Model object from mibitrans.transport.
        time (float): Point of time for the plot. By default, last point in time is plotted.
        **kwargs : Arguments to be passed to plt.plot_surface().

    Returns:
        ax (matplotlib.axes._axes.Axes) : Returns matplotlib axes object of plume plot.
    """
    _check_model_type(model, Domenico)
    t_pos = _time_check(model, time)

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(model.xxx[t_pos, :, :], model.yyy[t_pos, :, :], model.cxyt[t_pos, :, :], **kwargs)

    ax.view_init(elev=30, azim=310)
    ax.set_xlabel("Distance from source (m)")
    ax.set_ylabel("Distance from plume center (m)")
    ax.set_zlabel(r"Concentration [$g/m^{3}$]")

    return ax
