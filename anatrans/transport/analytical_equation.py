"""Author: Jorrit Bakker.

Module calculating the solution to the Domenico (1987) analytical model for different scenarios.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf
from scipy.special import erfc
from anatrans.data.check_input import CheckInput


class Transport:
    """Calculate analytical solution from variations of the Domenico analytical model."""
    def __init__(self,
                 parameters : dict,
                 dx : float = None,
                 dy : float = None,
                 dt : float = None,
                 verbose : bool = True
                 ) -> None:
        """Initialise the class, set up dimensions and output array."""
        self.pars = parameters
        self.verbose = verbose
        self.keys = self.pars.keys()

        self.ck = CheckInput(self.pars)
        par_success = self.ck.check_parameter()
        if not par_success:
            raise ValueError("Not all required input parameters are given.")

        val_success = self.ck.check_values()
        if not val_success:
            raise ValueError("Not all required input values are of the expected value or data type.")

        self.R = calculate_retardation(self.pars)
        self.v = calculate_flow_velocity(self.pars) / self.R
        self.ax, self.ay, self.az = calculate_dispersivity(self.pars)

        self.source_y = self.pars["c_source"][:,0]
        self.source_c = self.pars["c_source"][:,1]
        self.width = np.max(self.source_y) * 2

        # If dx, dy or dt are not given, determine a proper value based on modeled area dimensions.
        # If modeled area dimensions are not given, determine acceptable value based on flow velocity and retardation.
        self.dx = dx
        self.dy = dy
        self.dt = dt

        self.x = np.arange(0, parameters["l_model"], self.dx)
        # Center y dimension on the centerline of the plume
        self.y = np.zeros(2 * len(self.source_y) - 1)
        self.y[0:len(self.source_y)] = -self.source_y[::-1]
        self.y[len(self.source_y) - 1:] = self.source_y

        self.c0 = np.zeros(2 * len(self.source_c) - 1)
        self.c0[0:len(self.source_c)] = self.source_c[::-1]
        self.c0[len(self.source_c) - 1:] = self.source_c

        #self.y = np.arange(0, self.width, self.dy) - (self.width / 2)
        self.t = np.arange(0.00001, parameters["t_model"] + 0.00001, self.dt)
        #self.t[0] = self.t[0] + 0.9
        print(self.t[0] + 0.9)
        self.xxx = np.tile(self.x, (len(self.t), len(self.y), 1))
        self.yyy = np.tile(self.y[:, None], (len(self.t), 1, len(self.x)))
        self.ccc0 = np.tile(self.c0[:, None], (len(self.t), 1, len(self.x)))
        self.ttt = np.tile(self.t[:, None, None], (1, len(self.y), len(self.x)))

        self.cxyt = np.zeros(self.xxx.shape)

    def no_decay(self):
        """Calculate the Domenico analytical model without decay."""
        inf_erfc_x = (self.xxx - self.v * self.ttt) / (2 * np.sqrt(self.ax * self.v * self.ttt))
        erfc_x = erfc(inf_erfc_x)

        erf_y = (erf((self.yyy + self.width / 2) / (2 * np.sqrt(self.ay * self.xxx)))
                 - erf((self.yyy - self.width / 2) / (2 * np.sqrt(self.ay * self.xxx))))

        erf_z = (erf(self.pars["d_source"] / (2 * np.sqrt(self.az * self.xxx)))
                 - erf(-self.pars["d_source"] / (2 * np.sqrt(self.ay * self.xxx))))


        cxty = 1/8 * self.ccc0 * erfc_x * erf_y * erf_z

        return cxty, self.x, self.y, self.t

    def visualize_initial(self):
        plt.plot(self.y, self.c0)

def calculate_retardation(pars):
    """Give retardation factor depending on input parameters."""
    if "R" in pars.keys():
        r = pars["R"]
    else:
        r = 1 + (pars["rho"] / pars["n"]) * pars["Koc"] * pars["foc"]
    return r

def calculate_flow_velocity(pars):
    """Give flow velocity depending on input parameters."""
    if "v" in pars.keys():
        v = pars["v"]
    else:
        v = (pars["k"] * pars["i"]) / (pars["n"])
    return v

def calculate_dispersivity(pars):
    """Give dispersivity in each direction depending on input parameters."""
    if "alpha_x" in pars.keys():
        ax = pars["alpha_x"]
        ay = pars["alpha_y"]
        az = pars["alpha_z"]
    else:
        # Implement conversion from plume length to alpha x, y and z
        ax = 0.83 * np.log10(pars["lp"])**2.414
        ay = ax / 10
        az = ax / 100

    if ax < 1e-10:
        ax = 1e-10
    if ay < 1e-10:
        ay = 1e-10
    if az < 1e-10:
        az = 1e-10

    return ax, ay, az


