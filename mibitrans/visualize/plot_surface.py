import matplotlib.pyplot as plt
import numpy as np


class Plume():
    """Multi-dimensional plotting of contaminant plume."""
    def __init__(self, cxyt, x, y, t):
        """Initialize parameters."""
        self.cxyt = cxyt
        self.x = x
        self.y = y
        self.t = t
        self.xxx = np.tile(self.x, (len(self.t), len(self.y), 1))
        self.yyy = np.tile(self.y[:, None], (len(self.t), 1, len(self.x)))

    def surface(self, time = None, **kwargs):
        """Plot contaminant plume as a 3D surface plot."""
        if time is not None:
            time_pos = np.argmin(abs(self.t - time))
        else:
            time_pos = self.t[-1]

        fig, ax = plt.subplots(subplot_kw={"projection" : "3d"})
        ax.plot_surface(self.xxx[time_pos,:,:],
                        self.yyy[time_pos,:,:],
                        self.cxyt[time_pos,:,:],
                        **kwargs
                        )

        ax.view_init(elev=30, azim=310)
        ax.set_xlabel("Distance from source (m)")
        ax.set_ylabel("Distance from plume center (m)")
        ax.set_zlabel(r"Concentration [$g/m^{-3}$]")

        return(ax)

    def flat(self, time = None, **kwargs):
        """Plot contaminant plume as a 2D surface plot."""
        if time is not None:
            time_pos = np.argmin(abs(self.t - time))
        else:
            time_pos = self.t[-1]

        plt.pcolormesh(self.x,
                       self.y,
                       self.cxyt[time_pos,:,:],
                       **kwargs
                       )

        plt.xlabel("Distance from source (m)")
        plt.ylabel("Distance from plume center (m)")
        plt.colorbar(label = "Concentration (mg/L)")
