"""Author: Jorrit Bakker.

Module plotting a 3D matrix of contaminant plume concentrations as a line.
"""

import matplotlib.pyplot as plt
import numpy as np


class Lineplot():
    """Line plotting of contaminant plume."""
    def __init__(self, cxyt, x, y, t):
        """Initialize parameters."""
        self.cxyt = cxyt
        self.x = x
        self.y = y
        self.t = t

    def centerline(self, time = None, y_pos = 0, **kwargs):
        """Plot center of contaminant plume as a line, at a specified time and, optionally, y position.

        Args:
            time : Point of time for the plot. By default, last point in time is plotted.
            y_pos : y-position across the plume for the plot. By default, the center of the plume at y=0 is plotted.
            **kwargs : Arguments to be passed to plt.plot().

        Returns a line plot of the input plume as object.
        """
        if time is not None:
            time_pos = np.argmin(abs(self.t - time))
        else:
            time_pos = self.t[-1]

        if y_pos is not None:
            y_pos = np.argmin(abs(self.y - y_pos))
        else:
            y_pos = np.argmin(abs(self.y - 0))

        plot_array = self.cxyt[time_pos,y_pos,:]

        plt.plot(self.x,
                 plot_array,
                 **kwargs)
        plt.ylim((0,np.max(plot_array) + 1/8*np.max(plot_array)))
        plt.xlabel("Distance from source [m]")
        plt.ylabel("Concentration [mg/L]")
        plt.grid(True)


    def transverse(self, time = None, x_pos = None):
        """Plot across the contaminant plume as a line, at a specified time and x position ."""
        print("placeholder")

    def breakthrough(self, x_pos = None, y_pos = None):
        """Plot breakthrough curve of contaminant plume at a specified x and y position."""
        print("placeholder")
