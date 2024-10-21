"""Author: Jorrit Bakker.

Module calculating the solution to the Domenico (1987) analytical model for different scenarios.
"""

import numpy as np
from scipy.special import erf, erfc
from anatrans.data.check_input import CheckInput


class Transport():
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

        if "R" in self.keys:
            self.R = self.pars["R"]
        else:
            self.R = 1 + (self.pars["rho"] / self.pars["n"]) * self.pars["Koc"] * self.pars["foc"]

        if "v" in self.keys:
            self.v = self.pars["v"] / self.R
        else:
            self.v = (self.pars["k"] * self.pars["i"]) / (self.pars["n"] * self.R)

        if "alpha_x" in self.keys:
            self.ax = self.pars["alpha_x"]
            self.ay = self.pars["alpha_y"]
            self.az = self.pars["alpha_z"]
        else:
            # Implement conversion from plume length to alpha x, y and z
            self.ax = self.pars["lp"]

        if self.az == 0:
            self.az = 1e-10

        self.width = np.sum(self.pars["c_source"][:,0])
        # For now, initial concentration is the mean of the concentration across the source.
        self.c0 = np.mean(self.pars["c_source"][:,1])

        # If dx, dy or dt are not given, determine a proper value based on modeled area dimensions.
        # If modeled area dimensions are not given, determine acceptable value based on flow velocity and retardation.
        self.dx = dx
        self.dy = dy
        self.dt = dt

        self.x = np.arange(0, parameters["l_model"], self.dx)
        self.y = np.arange(0, self.width, self.dy) - (self.width / 2)
        self.t = np.arange(0, parameters["t_model"], self.dt)

        self.xxx = np.tile(self.x, (len(self.t), len(self.y), 1))
        self.yyy = np.tile(self.y[:, None], (len(self.t), 1, len(self.x)))
        self.ttt = np.tile(self.t[:, None, None], (1, len(self.y), len(self.x)))

        self.cxyt = np.zeros(self.xxx.shape)

    def no_decay(self):
        """Calculate the Domenico analytical model without decay."""
        erfc_x = erfc((self.xxx - self.v * self.ttt) / (2 * self.ax * self.v * self.ttt))
        #print(erfc_x)
        erf_y = (erf((self.yyy + self.width / 2) / (2 * np.sqrt(self.ay * self.xxx)))
                 - erf((self.yyy - self.width / 2) / (2 * np.sqrt(self.ay * self.xxx))))
        #print(erf_y)
        erf_z = (erf(self.pars["d_source"] / (2 * np.sqrt(self.az * self.xxx)))
                 - erf(-self.pars["d_source"] / (2 * np.sqrt(self.ay * self.xxx))))
        #print(erf_z)

        self.cxty = 1/8 * self.c0 * erfc_x * erf_y * erf_z

        return(self.cxty, self.x, self.y, self.t)
