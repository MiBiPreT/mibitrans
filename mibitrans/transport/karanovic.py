import numpy as np
from scipy.integrate import quad_vec
from mibitrans.transport.model_parent import Karanovic


class NoDecay(Karanovic):
    def __init__(
        self, hydrological_parameters, adsorption_parameters, source_parameters, model_parameters, verbose=False
    ):
        """Initialize object and run model.

        Args:
            hydrological_parameters (mibitrans.data.read.HydrologicalParameters) : Dataclass object containing
            hydrological parameters from HydrologicalParameters.
            adsorption_parameters (mibitrans.data.read.AdsorptionParameters) : Dataclass object containing adsorption
                parameters from AdsorptionParameters.
            source_parameters (mibitrans.data.read.SourceParameters) : Dataclass object containing source parameters
                from SourceParameters.
            model_parameters (mibitrans.data.read.ModelParameters) : Dataclass object containing model parameters from
                ModelParameters.
            verbose (bool, optional): Verbose mode. Defaults to False.

        Attributes:
            cxyt (np.ndarray) : Output array containing concentrations in model domain, in [g/m^3]. Indexed as [t,y,x]
            x (np.ndarray) : Discretized model x-dimension, in [m].
            y (np.ndarray) : Discretized model y-dimension, in [y].
            t (np.ndarray) : Discretized model t-dimension, in [days].
            c_source (np.ndarray) : Nett source zone concentrations, accounting for source superposition, in [g/m^3].
            vr (float) : Retarded groundwater flow velocity, in [m/d].
            k_source (float) : Source zone decay rate, in [1/days]

        Methods:
            sample : Give concentration at any given position and point in time, closest as discretization allows.

        Raises:
            TypeError : If input is not of the correct Dataclass.

        """
        super().__init__(hydrological_parameters, adsorption_parameters, source_parameters, model_parameters, verbose)
        self.integral_term = np.zeros(self.ttt.shape)
        self.error_size = np.zeros(len(self.t))
        self._calculate()

    def _calculate(self):
        with np.errstate(divide="ignore", invalid="ignore"):
            for i, time in enumerate(self.t):
                if self.verbose:
                    print("integrating for t =", time, "days")
                self.integral_term, self.error_size[i] = quad_vec(self._eq_integrand, 0, time)
            source_term = self._eq_source_term()
            self.cxyt[:, :, 1:] += self.integral_term * source_term
            self.cxyt[:, :, 0] = self._eq_source_zero()

    # def point(self, t, y, x):
    #     def integrand(t, y, x):
    #         x_part = np.exp(self.k_source * t - (x - self.rv * t) ** 2 / (4 * self.disp_x * t))
    #         div_term = 2 * np.sqrt(self.disp_y * t)
    #         y_part = erfc((y - self.y_source[0]) / div_term) - erfc((y + self.y_source[0]) / div_term)
    #         inner_term = self.src_pars.depth / (2 * np.sqrt(self.disp_z * t))
    #         z_part = erfc(-inner_term) - erfc(inner_term)
    #         term = 1 / (t ** (3 / 2)) * x_part * y_part * z_part
    #         return term
    #
    #     inter, _ = quad_vec(integrand, 0, t, args=(y, x))
    #     print(inter)
    #     src = self.c_source[0] * x / (8 * np.sqrt(np.pi * self.disp_x)) * np.exp(-self.k_source * self.t)
    #     print(src)
    #     c = inter * src
    #     return c
