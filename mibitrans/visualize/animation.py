import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation


def animate_1d(
    x_axis_parameter: np.ndarray, y_axis_parameter, time_parameter, y_names=None, y_colors=None, linestyle=None
):
    """Animate any type and number of input 1D plottables."""
    if not isinstance(y_axis_parameter, list):
        y_axis_parameter = [y_axis_parameter]
    fig, ax = plt.subplots()
    plot_bin = []
    max_conc = 0
    for i, y in enumerate(y_axis_parameter):
        line = ax.plot(x_axis_parameter, y[0, :], color=y_colors[i], label=y_names[i], linestyle=linestyle[i])[0]
        if np.max(y) > max_conc:
            max_conc = np.max(y)
        plot_bin.append(line)
    ax.set_ylim(bottom=0, top=max_conc + max_conc / 10)
    ax.set_xlabel("Distance from source [m]")
    ax.set_ylabel("Concentration [g/m3]")
    ax.legend()
    n_frames = len(time_parameter)

    def update(frame):
        """Update plot with values for the next time step in the animation."""
        for i, y in enumerate(y_axis_parameter):
            plot_bin[i].set_xdata(x_axis_parameter)
            plot_bin[i].set_ydata(y[frame, :])
        ax.set_title(f"Concentration distribution at t={time_parameter[frame]} days")
        return plot_bin

    ani = animation.FuncAnimation(fig=fig, func=update, frames=n_frames)
    return ani
