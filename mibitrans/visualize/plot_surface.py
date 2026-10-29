import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib import animation
import mibitrans
from mibitrans.data.check_input import check_model_type
from mibitrans.data.check_input import check_time_in_domain
from mibitrans.visualize.plot_line import _plot_title_generator
from mibitrans.visualize.plot_line import allowed_model_types

relative_conc_zlabel = r"Relative concentration ($C/C_0$)"
absolute_conc_zlabel = r"Concentration [g/$m^{3}$]"


def plume_core_fringe(
    model: mibitrans.transport.model_parent.Results,
    time: int | float = None,
    relative_concentration: bool = False,
    **kwargs,
):
    """Plot concentration distribution electron donor and electron acceptor together in a single contour plot.

    Args:
        model : Model Results object from mibitrans.transport, ran with core-fringe degradation.
        time (float): Point of time for the plot. Will show the closest time step to given value.
            By default, last point in time is plotted.
        relative_concentration (bool, optional) : If set to True, will plot concentrations relative to maximum source
            zone concentrations at t=0. By default, absolute concentrations are shown.
        **kwargs : Arguments for optional further plot customization:
            ed_nlevels (int) : Amount of contour levels for electron donor.
            ea_nlevels (int) : Amount of contour levels for electron acceptor.
            ed_cmap (str) : Colormap to use for electron donor contours.
            ea_cmap (str) : Colormap to use for electron acceptor contours.
            figsize (tuple) : Figure size of plot.
    """
    check_model_type(model, allowed_model_types())
    if model.mode != "core-fringe":
        raise ValueError(
            "Combined visualization for electron donor and electron acceptors is only possible for core-fringe models."
        )
    t_pos = check_time_in_domain(model, time)
    x = model.x
    y = model.y
    if relative_concentration:
        ed_cxyt, ea_cxyt = model.relative_cxyt[0], model.relative_cxyt[1]
    else:
        ed_cxyt, ea_cxyt = model.cxyt[0], model.cxyt[1]

    combi_plume = np.where(ed_cxyt <= 0, ed_cxyt - ea_cxyt, ed_cxyt)

    # MaxNLocator determines 'neat' intervals for colorbar automatically.
    tick_determination = ticker.MaxNLocator(nbins=12, steps=[1, 2, 4, 5, 10], min_n_ticks=6)
    # Can set amount of levels manually by specifying in kwargs
    if kwargs.get("ed_nlevels", None):
        ed_levels = np.linspace(0, np.max(ed_cxyt), kwargs["ed_nlevels"])
    else:
        ed_levels = np.array(tick_determination.tick_values(0, np.max(ed_cxyt)))
    if kwargs.get("ea_nlevels", None):
        ea_levels = np.linspace(-np.max(ea_cxyt), 0, kwargs["ea_nlevels"])
    else:
        ea_levels = np.array(tick_determination.tick_values(-np.max(ea_cxyt), 0))

    # Use default colormaps if colormaps not specified in kwargs
    ed_cmap = kwargs.get("ed_cmap", "YlOrBr")
    ea_cmap = kwargs.get("ea_cmap", "Greens_r")

    # Will use figsize of (16,9) if not specified differently in kwargs
    plt.figure(figsize=kwargs.get("figsize", (16, 9)), dpi=300, layout="constrained")
    plt.xlabel(r"Distance from source [$m$]")
    plt.ylabel(r"Distance from plume center [$m$]")
    plt.yticks([int(y[0]), 0, np.ceil(y[-1])], [np.round(y[0]), 0, np.round(y[-1])])

    ed_plume = plt.contourf(x, y, combi_plume[t_pos, :, :], levels=ed_levels, cmap=ed_cmap)
    ea_plume = plt.contourf(x, y, combi_plume[t_pos, :, :], levels=ea_levels, cmap=ea_cmap)
    # Line for boundary between ED and EA concentrations
    plume_max = plt.contour(x, y, combi_plume[t_pos, :, :], levels=[0], linewidths=2, colors="k")

    cbar_ed = plt.colorbar(ed_plume, label=r"Electron donor concentrations [$g/m^3$]", location="bottom", aspect=75)
    cbar_ea = plt.colorbar(ea_plume, label=r"Electron acceptor concentrations [$g/m^3$]", location="bottom", aspect=75)
    # Ensure that colorbar values for EA are positive
    cbar_ea.set_ticklabels([f"{round(abs(level), 5)}" for level in ea_plume.levels])

    bar_height = 0.01
    bar_width = 0.8
    bar_x = 0.1
    cbar_ea.ax.set_position([bar_x, 0, bar_width, bar_height])
    cbar_ed.ax.set_position([bar_x, 1, bar_width, bar_height])


def plume_2d(
    model: mibitrans.transport.model_parent.Results,
    time: int | float = None,
    relative_concentration: bool = False,
    animate: bool = False,
    plot_index: int | list[int] | None = None,
    **kwargs,
):
    """Plot contaminant plume as a 2D colormesh, at a specified time.

    Args:
        model : Model Results object from mibitrans.transport.
        time (float): Point of time for the plot. Will show the closest time step to given value.
            By default, last point in time is plotted.
        relative_concentration (bool, optional) : If set to True, will plot concentrations relative to maximum source
            zone concentrations at t=0. By default, absolute concentrations are shown.
        animate (bool, optional): If True, animation of contaminant plume until given time is shown. Default is
            False.
        plot_index (int, optional): Which concentration distribution to plot, if model has a list of multiple cxyt.
            As a zero-base index in the same order as decay rates were provided for chain-decay. For core-fringe,
            0-index contains electron donor distribution, 1-index contains electron acceptor distribution. Default is 0.
        **kwargs : Arguments to be passed to plt.pcolormesh().

    Returns a matrix plot of the input plume as object.
    """
    check_model_type(model, allowed_model_types())
    t_pos = check_time_in_domain(model, time)
    if isinstance(model.cxyt, list):
        model_concentration = model.cxyt[plot_index]
        z_label = absolute_conc_zlabel
    else:
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


def plume_3d(
    model: mibitrans.transport.model_parent.Results,
    time: int | float = None,
    relative_concentration: bool = False,
    animate: bool = False,
    plot_index: int | list[int] | None = None,
    **kwargs,
):
    """Plot contaminant plume as a 3D surface, at a specified time.

    Args:
        model : Model Results object from mibitrans.transport.
        time (float): Point of time for the plot. Will show the closest time step to given value.
            By default, last point in time is plotted.
        relative_concentration (bool, optional) : If set to True, will plot concentrations relative to maximum source
            zone concentrations at t=0. By default, absolute concentrations are shown.
        animate (bool, optional): If True, animation of contaminant plume until given time is shown. Default is
            False.
        plot_index (int, optional): Which concentration distribution to plot, if model has a list of multiple cxyt.
            As a zero-base index in the same order as decay rates were provided for chain-decay. For core-fringe,
            0-index contains electron donor distribution, 1-index contains electron acceptor distribution. Default is 0.
        **kwargs : Arguments to be passed to plt.plot_surface().

    Returns:
        ax (matplotlib.axes._axes.Axes) : Returns matplotlib axes object of plume plot.
    """
    check_model_type(model, allowed_model_types())
    t_pos = check_time_in_domain(model, time)
    if isinstance(model.cxyt, list):
        model_concentration = model.cxyt[plot_index]
        z_label = absolute_conc_zlabel
    else:
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
