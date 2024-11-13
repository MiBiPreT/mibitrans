"""Author: Jorrit Bakker.

Module calculating the solution to the Domenico (1987) analytical model for different scenarios.
"""

import numpy as np
from scipy.special import erf
from scipy.special import erfc
from anatrans.analysis.parameter_calculations import calculate_dispersivity
from anatrans.analysis.parameter_calculations import calculate_flow_velocity
from anatrans.analysis.parameter_calculations import calculate_linear_decay
from anatrans.analysis.parameter_calculations import calculate_retardation
from anatrans.analysis.parameter_calculations import calculate_source_decay
from anatrans.data.check_input import CheckInput


class Transport:
    """Calculate analytical solution from variations of the Domenico analytical model."""
    def __init__(self,
                 parameters : dict,
                 mode : str = None,
                 dx : float = None,
                 dy : float = None,
                 dt : float = None,
                 verbose : bool = True
                 ) -> None:
        """Initialise the class, set up dimensions and output array."""
        self.pars = parameters
        self.verbose = verbose
        self.mode = mode

        # Ensure that every parameter required for calculations is present in pars
        self.ck = CheckInput(self.pars, self.mode)
        par_success = self.ck.check_parameter()
        if not par_success:
            raise ValueError("Not all required input parameters are given.")

        # Ensure that every parameter required for calculations is of the correct value or data type.
        val_success = self.ck.check_values()
        if not val_success:
            raise ValueError("Not all required input values are of the expected value or data type.")

        self.R = calculate_retardation(self.pars)
        v = calculate_flow_velocity(self.pars)
        self.k_source = calculate_source_decay(self.pars)
        self.ax, self.ay, self.az = calculate_dispersivity(self.pars)

        if mode == "linear_decay":
            self.mu = calculate_linear_decay(self.pars)

        # For the rest of the calculation, retarded flow velocity is used.
        self.v = v / self.R

        # If dx, dy or dt are not given, determine a proper value based on modeled area dimensions.
        # If modeled area dimensions are not given, determine acceptable value based on flow velocity and retardation.
        self.dx = dx
        self.dy = dy
        self.dt = dt

        self.source_y = self.pars["c_source"][:,0]
        self.source_c = self.pars["c_source"][:,1]

        # Source zone concentration for superposition algorithm
        self.c0 = self.source_c[:-1] - self.source_c[1:]

        # Initialize space and time arrays
        self.x = np.arange(0, parameters["l_model"], self.dx)
        self.y = np.arange(-self.source_y[-1], self.source_y[-1] + self.dy, self.dy)
        self.t = np.arange(self.dt, parameters["t_model"] + self.dt, self.dt)

        # To allow array calculations, set space and time as a 3-dimensional array.
        self.xxx = np.tile(self.x, (len(self.t), len(self.y), 1))
        self.yyy = np.tile(self.y[:, None], (len(self.t), 1, len(self.x)))
        self.ttt = np.tile(self.t[:, None, None], (1, len(self.y), len(self.x)))
        self.cxyt = np.zeros(self.xxx.shape)

    def domenico(self):
        """Calculate the Domenico analytical model."""
        # terms for linear decay are calculated, terms are 1 for no decay and instant_reaction modes.
        if self.mode == "linear_decay":
            decay_sqrt = np.sqrt(1 + 4 * self.mu * self.ax / self.v)
            decay_exp = np.exp(self.xxx * (1 - decay_sqrt) / (self.ax * 2))
        elif self.mode == "instant_reaction":
            decay_sqrt = 1
            decay_exp = 1
        else:
            decay_sqrt = 1
            decay_exp = 1

        # Advection term and dispersion in the x direction
        erfc_x = erfc((self.xxx - self.v * self.ttt * decay_sqrt) / (2 * np.sqrt(self.ax * self.v * self.ttt)))

        # Dispersion term in the z direction
        erf_z = (erf(self.pars["d_source"] / (2 * np.sqrt(self.az * self.xxx)))
                 - erf(-self.pars["d_source"] / (2 * np.sqrt(self.az * self.xxx))))

        # Source decay term
        exp_source = np.exp(-self.k_source * (self.ttt - self.xxx / self.v))
        # Term can be max 1; can not have 'generation' of solute ahead of advection
        exp_source = np.where(exp_source > 1, 1, exp_source)

        # Superposition algorithm; calculate plume for each source zone separately, then add together.
        ccc0_source_list = [0] * len(self.c0)
        for i in range(len(self.c0)):
            # Dispersion term in the y direction
            erf_y = (erf((self.yyy + self.source_y[i+1]) / (2 * np.sqrt(self.ay * self.xxx)))
                     - erf((self.yyy - self.source_y[i+1]) / (2 * np.sqrt(self.ay * self.xxx))))

            ccc0_source_list[i] = self.c0[i] * exp_source
            cxyt = 1 / 8 * ccc0_source_list[i] * decay_exp * erfc_x * erf_y * erf_z
            self.cxyt += cxyt

        return self.cxyt, self.x, self.y, self.t
