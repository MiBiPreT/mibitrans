import warnings
import numpy as np
from scipy.special import erf
from scipy.special import erfc
import mibitrans.data.read
from mibitrans.data.check_input import validate_input_values


class Transport3D:
    """Parent class that for all 3-dimensional analytical solutions."""

    def __init__(
        self, hydrological_parameters, adsorption_parameters, source_parameters, model_parameters, verbose=False
    ):
        """Initialize parent class object.

        Args:
            hydrological_parameters (mibitrans.data.read.HydrologicalParameters) : Dataclass object containing
                hydrological parameters from HydrologicalParameters.
            adsorption_parameters (mibitrans.data.read.AdsorptionParameters) : Dataclass object containing
                adsorption parameters from AdsorptionParameters.
            source_parameters (mibitrans.data.read.SourceParameters) : Dataclass object containing source parameters
                from SourceParameters.
            model_parameters (mibitrans.data.read.ModelParameters) : Dataclass object containing model parameters from
                ModelParameters.
            verbose (bool, optional): Verbose mode. Defaults to False.
        """
        self.hyd_pars = hydrological_parameters
        self.ads_pars = adsorption_parameters
        self.src_pars = source_parameters
        self.mod_pars = model_parameters
        self.verbose = verbose

        # Check if input arguments are of the correct dataclass
        for key in self.__dict__.keys():
            if key != "verbose":
                self._check_input_dataclasses(key)

        # One-dimensional model domain arrays
        self.x = np.arange(0, self.mod_pars.model_length + self.mod_pars.dx, self.mod_pars.dx)
        self.y = self._calculate_y()
        self.t = np.arange(self.mod_pars.dt, self.mod_pars.model_time + self.mod_pars.dt, self.mod_pars.dt)

        # Three-dimensional model domain arrays
        self.xxx = np.tile(self.x, (len(self.t), len(self.y), 1))
        self.yyy = np.tile(self.y[:, None], (len(self.t), 1, len(self.x)))
        self.ttt = np.tile(self.t[:, None, None], (1, len(self.y), len(self.x)))
        # cxyt is concentration output array
        self.cxyt = np.zeros(self.xxx.shape)

        # Calculate retardation if not already specified in adsorption_parameters
        if self.ads_pars.retardation is None:
            self.ads_pars.calculate_retardation(self.hyd_pars.porosity)
        # Calculate retarded velocity
        self.rv = self.hyd_pars.velocity / self.ads_pars.retardation

        self.k_source = self._calculate_source_decay()
        self.y_source = self.src_pars.source_zone_boundary
        # Subtract outer source zones from inner source zones
        self.c_source = self.src_pars.source_zone_concentration.copy()
        self.c_source[:-1] = self.c_source[:-1] - self.c_source[1:]

    def sample(self, x_position, y_position, time):
        """Give concentration at any given position and point in time, closest as discretization allows.

        Args:
            x_position (float): x position in domain extent [m].
            y_position (float): y position in domain extent [m].
            time (float): time for which concentration is sampled [days].
            print_exact_location (bool, optional): If set to True, will print out exact location for which the
                concentration was determined. Defaults to False.

        Returns:
            concentration (float): concentration at given position and point in time [g/m^3]. Note that due to
                discretization, exact point of sampling can be up to half of a step size off in each dimension.

        """
        for par, value in locals().items():
            if par != "self":
                validate_input_values(par, value)

        # Save original 3d model domain and concentration arrays, so it can be restored afterward
        # Needed because _calculate uses xxx, yyy and ttt and would overwrite cxyt.
        save_x = self.xxx.copy()
        save_y = self.yyy.copy()
        save_t = self.ttt.copy()
        save_c = self.cxyt.copy()
        if hasattr(self, "cxyt_noBC"):
            save_c_noBC = self.cxyt_noBC.copy()
        self.xxx = np.array([x_position])
        self.yyy = np.array([y_position])
        self.ttt = np.array([time])
        self.cxyt = np.array([0.0])
        self._calculate()
        concentration = self.cxyt[0]
        self.xxx = save_x
        self.yyy = save_y
        self.ttt = save_t
        self.cxyt = save_c
        if hasattr(self, "cxyt_noBC"):
            self.cxyt_noBC = save_c_noBC
        return concentration

    def _calculate_source_decay(self, biodegradation_capacity=0):
        """Calculate source decay, for instant_reaction, biodegradation_capacity is required."""
        if isinstance(self.src_pars.total_mass, (float, int)):
            y_src = np.zeros(len(self.src_pars.source_zone_boundary) + 1)
            y_src[1:] = self.src_pars.source_zone_boundary
            c_src = self.src_pars.source_zone_concentration
            Q = self.hyd_pars.velocity * self.hyd_pars.porosity * self.src_pars.depth * np.max(y_src) * 2

            weighted_conc = np.zeros(len(self.src_pars.source_zone_boundary))
            for i in range(len(self.src_pars.source_zone_boundary)):
                weighted_conc[i] = (y_src[i + 1] - y_src[i]) * c_src[i]

            c0_avg = biodegradation_capacity + np.sum(weighted_conc) / np.max(y_src)
            k_source = Q * c0_avg / self.src_pars.total_mass
        # If source mass is not a float, it is an infinite source, therefore, no source decay takes place.
        else:
            k_source = 0

        return k_source

    def _check_input_dataclasses(self, expected_class):
        """Check if input parameters are the correct dataclasses. Raise an error if not."""
        dataclass_dict = {
            "hyd_pars": mibitrans.data.read.HydrologicalParameters,
            "ads_pars": mibitrans.data.read.AdsorptionParameters,
            "deg_pars": mibitrans.data.read.DegradationParameters,
            "src_pars": mibitrans.data.read.SourceParameters,
            "mod_pars": mibitrans.data.read.ModelParameters,
        }

        rename_dict = {
            "hyd_pars": "hydrological_parameters",
            "ads_pars": "adsorption_parameters",
            "deg_pars": "degradation_parameters",
            "src_pars": "source_parameters",
            "mod_pars": "model_parameters",
        }

        if not isinstance(self.__dict__[expected_class], dataclass_dict[expected_class]):
            raise TypeError(
                f"Input argument {rename_dict[expected_class]} should be {dataclass_dict[expected_class]}, "
                f"but is {type(self.__dict__[expected_class])} instead."
            )

    def _calculate_y(self):
        """Calculate y-direction discretization."""
        if self.mod_pars.model_width >= 2 * self.src_pars.source_zone_boundary[-1]:
            y = np.arange(
                -self.mod_pars.model_width / 2, self.mod_pars.model_width / 2 + self.mod_pars.dy, self.mod_pars.dy
            )
        else:
            y = np.arange(
                -self.src_pars.source_zone_boundary[-1],
                self.src_pars.source_zone_boundary[-1] + self.mod_pars.dy,
                self.mod_pars.dy,
            )
            warnings.warn(
                "Source zone boundary is larger than model width. Model width adjusted to fit entire source zone."
            )
        return y


class Domenico(Transport3D):
    """Parent class that for all analytical solutions using on the Domenico (1987) analytical model.

    Domenico, P. A. (1987). An analytical model for multidimensional transport of a decaying contaminant species.
    Journal of Hydrology, 91(1-2), 49-58.
    """

    def _eq_x_term(self, decay_sqrt=1):
        return erfc(
            (self.xxx - self.hyd_pars.velocity * self.ttt * decay_sqrt)
            / (2 * np.sqrt(self.hyd_pars.alpha_x * self.hyd_pars.velocity * self.ttt))
        )

    def _eq_additional_x(self):
        return np.exp(self.xxx * self.hyd_pars.velocity / (self.hyd_pars.alpha_x * self.hyd_pars.velocity)) * (
            erfc(
                self.xxx
                + self.hyd_pars.velocity
                * self.ttt
                / (2 * np.sqrt(self.hyd_pars.alpha_x * self.hyd_pars.velocity * self.ttt))
            )
        )

    def _eq_z_term(self):
        inner_term = self.src_pars.depth / (2 * np.sqrt(self.hyd_pars.alpha_z * self.xxx))
        return erf(inner_term) - erf(-inner_term)

    def _eq_source_decay(self):
        term = np.exp(-self.k_source * (self.ttt - self.xxx / self.hyd_pars.velocity))
        # Term can be max 1; can not have 'generation' of solute ahead of advection.
        return np.where(term > 1, 1, term)

    def _eq_y_term(self, i):
        div_term = 2 * np.sqrt(self.hyd_pars.alpha_y * self.xxx)
        term = erf((self.yyy + self.y_source[i]) / div_term) - erf((self.yyy - self.y_source[i]) / div_term)
        term[np.isnan(term)] = 0
        return term


class Karanovic(Transport3D):
    def __init__(
        self,
        hydrological_parameters,
        adsorption_parameters,
        source_parameters,
        model_parameters,
        verbose=False,
    ):
        """Initialize object and run model.

        Args:
            hydrological_parameters (mibitrans.data.read.HydrologicalParameters) : Dataclass object containing
                hydrological parameters from HydrologicalParameters.
            adsorption_parameters (mibitrans.data.read.AdsorptionParameters) : Dataclass object containing adsorption
                parameters from AdsorptionParameters.
            degradation_parameters (mibitrans.data.read.DegradationParameters) : Dataclass object containing degradation
                parameters from DegradationParameters.
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
        self.disp_x = self.hyd_pars.alpha_x * self.rv  # Original equation also considers diffusion
        self.disp_y = self.hyd_pars.alpha_y * self.rv
        self.disp_z = self.hyd_pars.alpha_z * self.rv

    def _eq_integrand(self, t):
        term = 1 / (t ** (3 / 2)) * self._eq_x_exp_term(t) * self._eq_y_term(t) * self._eq_z_term(t)
        term[np.isnan(term)] = 0
        return term

    def _eq_x_exp_term(self, t):
        term = np.exp(self.k_source * t - (self.xxx[0, :, 1:] - self.rv * t) ** 2 / (4 * self.disp_x * t))
        term[np.isnan(term)] = 0
        return term

    def _eq_y_term(self, t):
        div_term = 2 * np.sqrt(self.disp_y * t)
        term = erfc((self.yyy[0, :, 1:] - self.y_source[0]) / div_term) - erfc(
            (self.yyy[0, :, 1:] + self.y_source[0]) / div_term
        )
        term[np.isnan(term)] = 0
        return term

    def _eq_z_term(self, t):
        inner_term = self.src_pars.depth / (2 * np.sqrt(self.disp_z * t))
        if inner_term == np.inf:
            inner_term = 2
        return erfc(-inner_term) - erfc(inner_term)

    def _eq_source_term(self):
        return (
            self.c_source[0]
            * self.xxx[:, :, 1:]
            / (8 * np.sqrt(np.pi * self.disp_x))
            * np.exp(-self.k_source * self.ttt[:, :, 1:])
        )

    def _eq_source_zero(self):
        return self.c_source[0] * np.exp(-self.k_source * self.ttt[:, :, 0])
