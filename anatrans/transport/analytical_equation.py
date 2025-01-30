"""Author: Jorrit Bakker.

Module calculating the solution to the Domenico (1987) analytical model for different scenarios.
"""

import numpy as np
from scipy.special import erf
from scipy.special import erfc
from anatrans.analysis.parameter_calculations import calculate_biodegradation_capacity
from anatrans.analysis.parameter_calculations import calculate_dispersivity
from anatrans.analysis.parameter_calculations import calculate_flow_velocity
from anatrans.analysis.parameter_calculations import calculate_linear_decay
from anatrans.analysis.parameter_calculations import calculate_retardation
from anatrans.analysis.parameter_calculations import calculate_source_decay
from anatrans.analysis.parameter_calculations import calculate_source_decay_instant
from anatrans.data.check_input import CheckInput


class Transport:
    """Calculate analytical solution from variations of the Domenico analytical model."""
    def __init__(self,
                 parameters : dict,
                 mode : str = "no_decay",
                 dx : float = 0,
                 dy : float = 0,
                 dt : float = 0,
                 verbose : bool = False
                 ) -> None:
        """Initialise the class, set up dimensions and output array.

        Args:
            parameters (dict) : Dictionary with transport parameters.
            mode (str) : Type of analytical model to be used. Default is no decay model.
            dx (float) : Model step size in x direction. If left empty,
            reasonable value will be calculated based on modeled area length. Default is None.
            dy (float) : Model step size in y direction. If left empty,
            reasonable value will be calculated based on modeled area width. Default is None.
            dt (float) : Model step size for time. If left empty,
            reasonable value will be calculated based on simulation time. Default is None.
            verbose (bool, optional): Verbose mode. Defaults to False.

        Raises:
            ValueError : If input parameters are incomplete or of incorrect datatype/value.

        """
        self.pars = parameters
        self.mode = mode
        self.verbose = verbose

        if verbose:
            if self.mode == "no_decay":
                print("Initializing parameters for the no decay model.")
            elif self.mode == "linear_decay":
                print("Initializing parameters for the linear decay model.")
            elif self.mode == "instant_reaction":
                print("Initializing parameters for the instantaneous reaction model.")
            else:
                print("Mode name is not recognized, initializing decay model as default. Valid model names are:",
                      " 'no_decay', 'linear_decay', 'instant_reaction'.")

        if verbose:
            print("Checking input parameters...")

        # Ensure that every parameter required for calculations is present in pars
        self.ck = CheckInput(self.pars, self.mode, verbose = self.verbose)
        par_success = self.ck.check_parameter()
        if not par_success:
            raise ValueError("Not all required input parameters are given.")

        # Ensure that every parameter required for calculations is of the correct value or data type.
        val_success = self.ck.check_values()
        if not val_success:
            raise ValueError("Not all required input values are of the expected value or data type.")

        if verbose:
            print("Calculating transport parameters...")

        self.R = calculate_retardation(self.pars)
        v = calculate_flow_velocity(self.pars)
        self.pars["v"] = v

        # No source decay occurs when source mass is infinite
        if self.pars["m_total"] in ["inf", "infinite"]:
            self.k_source = 0
        else:
            self.k_source = calculate_source_decay(self.pars)
        self.ax, self.ay, self.az = calculate_dispersivity(self.pars)

        if mode == "linear_decay":
            self.mu = calculate_linear_decay(self.pars)
        else:
            self.mu = 0

        # For the rest of the calculation, retarded flow velocity is used.
        self.v = v / self.R

        self.dx = dx
        self.dy = dy
        self.dt = dt

        # If dx, dy or dt are not given, determine a proper value based on modeled area dimensions.
        if not isinstance(dx, (float, int)) or (dx <= 0.0):
            self.dx = self.pars["l_model"] / 100
            if verbose:
                print("dx was not given or of wrong type/value, taken as fraction of model length instead.")
        if not isinstance(dy, (float, int)) or (dy <= 0.0):
            self.dy = self.pars["w_model"] / 50
            if verbose:
                print("dy was not given or of wrong type/value, taken as fraction of model width instead.")
        if not isinstance(dt, (float, int)) or (dt <= 0.0):
            self.dt = self.pars["t_model"] / 10
            if verbose:
                print("dt was not given or of wrong type/value, taken as fraction of model time instead.")

        self.source_y = self.pars["c_source"][:,0]
        self.source_c = self.pars["c_source"][:,1]

        # Source zone concentration for superposition algorithm
        self.c0 = self.source_c[:-1] - self.source_c[1:]

        # Calculate and add biodegradation capacity to the outer plume
        if mode == "instant_reaction":
            self.biodeg_cap = calculate_biodegradation_capacity(self.pars)
            self.c0[-1] += self.biodeg_cap

            # No source decay occurs when source mass is infinite
            if self.pars["m_total"] in ["inf", "infinite"]:
                self.k_source_instant = 0
            else:
                self.k_source_instant = calculate_source_decay_instant(self.pars, self.biodeg_cap)

        else:
            self.biodeg_cap = 0

        if verbose:
            print("Setting up dimensional grids...")

        # Initialize space and time arrays
        self.x = np.arange(0, parameters["l_model"] + self.dx, self.dx)

        if parameters["w_model"] > 2 * self.source_y[-1]:
            self.y = np.arange(-parameters["w_model"] / 2, parameters["w_model"] / 2, self.dy)
        else:
            self.y = np.arange(-self.source_y[-1], self.source_y[-1] + self.dy, self.dy)

        self.t = np.arange(self.dt, parameters["t_model"] + self.dt, self.dt)

        # To allow array calculations, set space and time as a 3-dimensional array.
        self.xxx = np.tile(self.x, (len(self.t), len(self.y), 1))
        self.yyy = np.tile(self.y[:, None], (len(self.t), 1, len(self.x)))
        self.ttt = np.tile(self.t[:, None, None], (1, len(self.y), len(self.x)))
        self.cxyt = np.zeros(self.xxx.shape)

        self.cxyt_noBC = 0

        if verbose:
            print("Done with initializing!")

    def domenico(self):
        """Calculate the Domenico analytical model.

        Returns:
            cxyt (np.ndarray) : 3-dimensional array with values representing contaminant concentrations. Outer dimension
            represents time, middle dimension represents y-direction and inner dimension represents x-direction.
            x (np.ndarray) : Array with values representing positions along the x-direction.
            y (np.ndarray) : Array with values representing positions along the y-direction.
            t (np.ndarray) : Array with values representing points in time.
        """
        # Ignore divide by zero warnings; they are a result of working with x and y positions of 0.
        with np.errstate(divide="ignore", invalid="ignore"):
            # terms for linear decay are calculated, terms are 1 for no decay and instant_reaction modes, as mu = 0.
            decay_sqrt = np.sqrt(1 + 4 * self.mu * self.ax / self.v)
            decay_term = np.exp(self.xxx * (1 - decay_sqrt) / (self.ax * 2))

            # Advection and dispersion term in the x direction
            x_term = erfc((self.xxx - self.v * self.ttt * decay_sqrt) / (2 * np.sqrt(self.ax * self.v * self.ttt)))

            # Additional advection and dispersion term in x direction for small times, when there is no linear decay
            if self.mode != "linear_decay":
                additional_x = (np.exp(self.xxx * self.v / (self.ax * self.v))
                                * erfc(self.xxx + self.v * self.ttt / (2 * np.sqrt(self.ax * self.v * self.ttt))))
                x_term += additional_x

            # Dispersion term in the z direction
            z_term = (erf(self.pars["d_source"] / (2 * np.sqrt(self.az * self.xxx)))
                     - erf(-self.pars["d_source"] / (2 * np.sqrt(self.az * self.xxx))))

            if self.mode == "instant_reaction":
                sourcedecay_term = np.exp(-self.k_source_instant * (self.ttt - self.xxx / self.v))
            else:
                # Source decay term
                sourcedecay_term = np.exp(-self.k_source * (self.ttt - self.xxx / self.v))
            # Term can be max 1; can not have 'generation' of solute ahead of advection
            sourcedecay_term = np.where(sourcedecay_term > 1, 1, sourcedecay_term)

            # Superposition algorithm; calculate plume for each source zone separately, then add together.
            ccc0_source_list = [0] * len(self.c0)
            for i in range(len(self.c0)):
                # Dispersion term in the y direction
                y_term = (erf((self.yyy + self.source_y[i+1]) / (2 * np.sqrt(self.ay * self.xxx)))
                         - erf((self.yyy - self.source_y[i+1]) / (2 * np.sqrt(self.ay * self.xxx))))

                y_term[np.isnan(y_term)] = 0

                # Correct source zone concentrations for source decay
                ccc0_source_list[i] = self.c0[i] * sourcedecay_term

                cxyt = 1 / 8 * ccc0_source_list[i] * decay_term * x_term * y_term * z_term

                self.cxyt += cxyt

        # Substract biodegradation capacity, this can cause low concentration areas to become < 0. Therefore, set
        # negative values to 0.
        self.cxyt_noBC = self.cxyt.copy()
        self.cxyt -= self.biodeg_cap
        self.cxyt = np.where(self.cxyt < 0, 0, self.cxyt)

        return self.cxyt, self.x, self.y, self.t
