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
        # with np.errstate(divide="ignore", invalid="ignore"):
        for sz in range(len(self.c_source)):
            if self.verbose:
                print("integrating for source zone ", sz)
            for j in range(len(self.t)):
                if self.verbose:
                    print("integrating for t =", self.t[j], "days")
                if j == 0:
                    lower_bound = 0
                else:
                    lower_bound = self.t[j - 1]
                upper_bound = self.t[j]
                self.integral_term[j, :, 1:], self.error_size[j] = quad_vec(
                    self._eq_integrand, lower_bound, upper_bound, limit=10000 / len(self.t), args=(sz,)
                )
            self.integral_sum = np.cumsum(self.integral_term, axis=0)
            source_term = self._eq_source_term(sz)
            self.cxyt[:, :, 1:] += self.integral_sum[:, :, 1:] * source_term
            self.cxyt[:, :, 0] += self._eq_source_zero(sz)
