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
                 initial : str = "direct",
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

        self.source_y = self.pars["c_source"][:,0]
        self.source_c = self.pars["c_source"][:,1]

        self.R = calculate_retardation(self.pars)
        self.v = calculate_flow_velocity(self.pars)

        self.source_decay()

        self.v = self.v / self.R
        self.ax, self.ay, self.az = calculate_dispersivity(self.pars)

        self.width = np.max(self.source_y) * 2

        # If dx, dy or dt are not given, determine a proper value based on modeled area dimensions.
        # If modeled area dimensions are not given, determine acceptable value based on flow velocity and retardation.
        self.dx = dx
        self.dy = dy
        self.dt = dt

        self.x = np.arange(0, parameters["l_model"], self.dx)

        if initial == "direct":
            self.initial_direct()
        elif initial == "superposition":
            self.initial_superposition()
        else:
            print("initial method is unknown, defaulting to direct")
            self.initial_direct()


        # Center y dimension on the centerline of the plume
        # self.y = np.zeros(2 * len(self.source_y) - 1)
        # self.y[0:len(self.source_y)] = -self.source_y[::-1]
        # self.y[len(self.source_y) - 1:] = self.source_y
        #
        # self.c0 = np.zeros(2 * len(self.source_c) - 1)
        # self.c0[0:len(self.source_c)] = self.source_c[::-1]
        # self.c0[len(self.source_c) - 1:] = self.source_c

        #self.y = np.arange(0, self.width, self.dy) - (self.width / 2)
        self.t = np.arange(0.00001, parameters["t_model"] + 0.00001 + self.dt, self.dt)
        #self.t[0] = self.t[0] + 0.9
        print(self.t[0] + 0.9)
        self.xxx = np.tile(self.x, (len(self.t), len(self.y), 1))
        self.yyy = np.tile(self.y[:, None], (len(self.t), 1, len(self.x)))

        self.ttt = np.tile(self.t[:, None, None], (1, len(self.y), len(self.x)))

        if self.c0.ndim == 2:
            self.ccc0_list = [0] * len(self.c0[0,:])
            print(len(self.ccc0_list))
            print(self.c0[:,0])
            for i in range(len(self.ccc0_list)):
                ccc0 = np.tile(self.c0[:, i][:, None], (len(self.t), 1, len(self.x)))
                self.ccc0_list[i] = ccc0
            print(ccc0.shape)
        elif self.c0.ndim == 1:
            self.width = np.tile(self.source_y)
            self.ccc0 = np.tile(self.c0[:, None], (len(self.t), 1, len(self.x)))
        else:
            print("something went very wrong")
        self.cxyt = np.zeros(self.xxx.shape)

    def no_decay(self):
        """Calculate the Domenico analytical model without decay."""

        inf_erfc_x = (self.xxx - self.v * self.ttt) / (2 * np.sqrt(self.ax * self.v * self.ttt))
        erfc_x = erfc(inf_erfc_x)

        erf_z = (erf(self.pars["d_source"] / (2 * np.sqrt(self.az * self.xxx)))
                 - erf(-self.pars["d_source"] / (2 * np.sqrt(self.az * self.xxx))))

        exp_source = np.exp(-self.k_source * (self.ttt - self.xxx / self.v))
        exp_source = np.where(exp_source > 1, 1, exp_source)

        if self.c0.ndim == 2:
            ccc0_source_list = [0] * len(self.ccc0_list)
            for i in range(len(self.c0[0,:])):
                erf_y = (erf((self.yyy + self.source_y[i+1]) / (2 * np.sqrt(self.ay * self.xxx)))
                         - erf((self.yyy - self.source_y[i+1]) / (2 * np.sqrt(self.ay * self.xxx))))

                ccc0_source_list[i] = self.ccc0_list[i] * exp_source
                cxyt = 1 / 8 * ccc0_source_list[i] * erfc_x * erf_y * erf_z
                self.cxyt += cxyt
            ccc0_array = np.asarray(ccc0_source_list)
        elif self.c0.ndim == 1:
            erf_y = (erf((self.yyy / 2) / (2 * np.sqrt(self.ay * self.xxx)))
                     - erf((self.yyy / 2) / (2 * np.sqrt(self.ay * self.xxx))))
            self.cxyt = 1 / 8 * self.ccc0 * exp_source * erfc_x * erf_y * erf_z

        else:
            print("something went very wrong")

        return self.cxyt, self.x, self.y, self.t, erfc_x

    def initial_superposition(self):
        """Function that calculates initial condition as a superposition of individual plumes varying in width
         with constant concentration.
        """
        self.y = np.arange(-self.source_y[-1], self.source_y[-1] + self.dy, self.dy)

        self.c0 = np.zeros((len(self.y), len(self.source_c) - 1))
        for i in range(len(self.source_c) - 1):
            con_value = self.source_c[i] - self.source_c[i+1]
            self.c0[:, i] = con_value

    def initial_direct(self):
        """Function that calculates initial condition as a plume with varying concentration across its width."""
        self.y = np.zeros(2 * len(self.source_y) - 1)
        self.y[0:len(self.source_y)] = -self.source_y[::-1]
        self.y[len(self.source_y) - 1:] = self.source_y

        self.c0 = np.zeros(2 * len(self.source_c) - 1)
        self.c0[0:len(self.source_c)] = self.source_c[::-1]
        self.c0[len(self.source_c) - 1:] = self.source_c

    def source_decay(self):
        """Function that calculates the source zone decay constant."""
        Q = self.v * self.pars["n"] * self.pars["d_source"] * np.max(self.source_y) * 2 * 1000

        C0_avg = 0
        for i in range(len(self.source_y) - 1):
            if i == 0:
                yc = self.source_y[i + 1] * self.source_c[i] * 2
            else:
                yc = (self.source_y[i+1] - self.source_y[i]) * self.source_c[i] * 2
            C0_avg += yc
        C0_avg = C0_avg / (np.max(self.source_y) * 2)

        self.k_source = Q * C0_avg / (self.pars["m_total"] * 1e6)


    def visualize_initial(self):
        """Temporary function to check if initial condition is correct."""
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


