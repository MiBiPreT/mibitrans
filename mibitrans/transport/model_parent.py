import copy
import warnings
from abc import ABC
from abc import abstractmethod
import numpy as np
from scipy.integrate import quad
from scipy.integrate import quad_vec
from scipy.special import erf
from scipy.special import erfc
import mibitrans.data.parameters
from mibitrans.data.check_input import validate_input_values
from mibitrans.visualize import plot_line as pline
from mibitrans.visualize import plot_surface as psurf


class Transport3D:
    """Parent class for all 3-dimensional analytical solutions."""

    def __init__(
        self, hydrological_parameters, attenuation_parameters, source_parameters, model_parameters, verbose=False
    ):
        """Initialize parent class object.

        Args:
            hydrological_parameters (mibitrans.data.parameters.HydrologicalParameters) : Dataclass object containing
                hydrological parameters from HydrologicalParameters.
            attenuation_parameters (mibitrans.data.read.AttenuationParameters) : Dataclass object containing adsorption,
                degradation and diffusion parameters from AttenuationParameters.
            source_parameters (mibitrans.data.read.SourceParameters) : Dataclass object containing source parameters
                from SourceParameters.
            model_parameters (mibitrans.data.read.ModelParameters) : Dataclass object containing model parameters from
                ModelParameters.
            verbose (bool, optional): Verbose mode. Defaults to False.
        """
        # Check if input arguments are of the correct dataclass
        for key, value in locals().items():
            if key not in ["self", "verbose"]:
                self._check_input_dataclasses(key, value)

        self._hyd_pars = copy.copy(hydrological_parameters)
        self._att_pars = copy.copy(attenuation_parameters)
        self._src_pars = copy.copy(source_parameters)
        self._mod_pars = copy.copy(model_parameters)

        self.verbose = verbose

        self._observe_input_dataclass_change()

        self.has_run = False
        self.initialized = False
        self._pre_run_initialization_parameters()

    @property
    def hydrological_parameters(self):
        """Rename to shorthand form of hydrological_parameters inside class for ease of use."""
        return self._hyd_pars

    @hydrological_parameters.setter
    def hydrological_parameters(self, value):
        self._hyd_pars = copy.copy(value)
        self._check_and_reset_when_input_dataclass_change("hydrological_parameters", value)

    @property
    def attenuation_parameters(self):
        """Rename to shorthand form of attenuation_parameters inside class for ease of use."""
        return self._att_pars

    @attenuation_parameters.setter
    def attenuation_parameters(self, value):
        self._att_pars = copy.copy(value)
        self._check_and_reset_when_input_dataclass_change("attenuation_parameters", value)

    @property
    def source_parameters(self):
        """Rename to shorthand form of source_parameters inside class for ease of use."""
        return self._src_pars

    @source_parameters.setter
    def source_parameters(self, value):
        self._src_pars = copy.copy(value)
        self._check_and_reset_when_input_dataclass_change("source_parameters", value)

    @property
    def model_parameters(self):
        """Rename to shorthand form of model_parameters inside class for ease of use."""
        return self._mod_pars

    @model_parameters.setter
    def model_parameters(self, value):
        self._mod_pars = copy.copy(value)
        self._check_and_reset_when_input_dataclass_change("model_parameters", value)

    def _observe_input_dataclass_change(self):
        self._hyd_pars._on_change = lambda: self._check_and_reset_when_input_dataclass_change(
            "hydrological_parameters", self._hyd_pars
        )
        self._att_pars._on_change = lambda: self._check_and_reset_when_input_dataclass_change(
            "attenuation_parameters", self._att_pars
        )
        self._src_pars._on_change = lambda: self._check_and_reset_when_input_dataclass_change(
            "source_parameters", self._src_pars
        )
        self._mod_pars._on_change = lambda: self._check_and_reset_when_input_dataclass_change(
            "model_parameters", self._mod_pars
        )

    def _check_and_reset_when_input_dataclass_change(self, key, value):
        self._check_input_dataclasses(key, value)
        self.initialized  = False
        if self.has_run:
            self.cxyt = np.zeros(self.xxx.shape)
            self.has_run = False
            if self.verbose:
                print(f"Parameter '{key}' has changed — resetting output")

    def _pre_run_initialization_parameters(self):
        # One-dimensional model domain arrays
        self.x = np.arange(0, self._mod_pars.model_length + self._mod_pars.dx, self._mod_pars.dx)
        self.y = self._calculate_y_discretization()
        self.t = np.arange(self._mod_pars.dt, self._mod_pars.model_time + self._mod_pars.dt, self._mod_pars.dt)

        # Three-dimensional model domain arrays
        self.xxx = np.tile(self.x, (len(self.t), len(self.y), 1))
        self.yyy = np.tile(self.y[:, None], (len(self.t), 1, len(self.x)))
        self.ttt = np.tile(self.t[:, None, None], (1, len(self.y), len(self.x)))

        self.rv = self._hyd_pars.velocity / self._att_pars.retardation

        # cxyt is concentration output array
        self.cxyt = np.zeros(self.xxx.shape)

        # Calculate retardation if not already specified in adsorption_parameters
        if (
            self._att_pars.bulk_density is not None
            and self._att_pars.partition_coefficient is not None
            and self._att_pars.fraction_organic_carbon is not None
        ):
            self._att_pars.calculate_retardation(self._hyd_pars.porosity)

        self.k_source = self._calculate_source_decay()
        self.y_source = self._src_pars.source_zone_boundary
        # Subtract outer source zones from inner source zones
        self.c_source = self._src_pars.source_zone_concentration.copy()
        self.c_source[:-1] = self.c_source[:-1] - self.c_source[1:]

        self.initialized = True

    def _calculate_source_decay(self, biodegradation_capacity=0):
        """Calculate source decay, for instant_reaction, biodegradation_capacity is required."""
        if isinstance(self._src_pars.total_mass, (float, int)):
            y_src = np.zeros(len(self._src_pars.source_zone_boundary) + 1)
            y_src[1:] = self._src_pars.source_zone_boundary
            c_src = self._src_pars.source_zone_concentration
            Q = self._hyd_pars.velocity * self._hyd_pars.porosity * self._src_pars.depth * np.max(y_src) * 2

            weighted_conc = np.zeros(len(self._src_pars.source_zone_boundary))
            for i in range(len(self._src_pars.source_zone_boundary)):
                weighted_conc[i] = (y_src[i + 1] - y_src[i]) * c_src[i]

            c0_avg = biodegradation_capacity + np.sum(weighted_conc) / np.max(y_src)
            k_source = Q * c0_avg / self._src_pars.total_mass
        # If source mass is not a float, it is an infinite source, therefore, no source decay takes place.
        else:
            k_source = 0

        return k_source

    def _check_input_dataclasses(self, key, value):
        """Check if input parameters are the correct dataclasses. Raise an error if not."""
        dataclass_dict = {
            "hydrological_parameters": mibitrans.data.parameters.HydrologicalParameters,
            "attenuation_parameters": mibitrans.data.parameters.AttenuationParameters,
            "source_parameters": mibitrans.data.parameters.SourceParameters,
            "model_parameters": mibitrans.data.parameters.ModelParameters,
        }

        if not isinstance(value, dataclass_dict[key]):
            raise TypeError(f"Input argument {key} should be {dataclass_dict[key]}, but is {type(value)} instead.")

    def _calculate_y_discretization(self):
        """Calculate y-direction discretization."""
        if self._mod_pars.model_width >= 2 * self._src_pars.source_zone_boundary[-1]:
            y = np.arange(
                -self._mod_pars.model_width / 2, self._mod_pars.model_width / 2 + self._mod_pars.dy, self._mod_pars.dy
            )
        else:
            y = np.arange(
                -self._src_pars.source_zone_boundary[-1],
                self._src_pars.source_zone_boundary[-1] + self._mod_pars.dy,
                self._mod_pars.dy,
            )
            warnings.warn(
                "Source zone boundary is larger than model width. Model width adjusted to fit entire source zone."
            )
        return y

    def sample(self, x_position, y_position, time):
        """Give concentration at any given position and point in time.

        Args:
            x_position (float): x position in domain extent [m].
            y_position (float): y position in domain extent [m].
            time (float): time for which concentration is sampled [days].

        Returns:
            concentration (float): concentration at given position and point in time [g/m^3].

        """
        for par, value in locals().items():
            if par != "self":
                validate_input_values(par, value)
        if not self.has_run and not self.initialized:
            self._pre_run_initialization_parameters()

        if hasattr(self, "cxyt_noBC"):
            save_c_noBC = self.cxyt_noBC.copy()
        x = np.array([x_position])
        y = np.array([y_position])
        t = np.array([time])
        concentration = self._calculate_cxyt(x, y, t)[0]
        if hasattr(self, "cxyt_noBC"):
            self.cxyt_noBC = save_c_noBC
        return concentration

    def centerline(self, y_position=0, time=None, legend_names=None, animate=False, **kwargs):
        """Plot center of contaminant plume of this model, at a specified time and y position.

        Args:
            y_position (float): y-position across the plume (transverse horizontal direction) for the plot.
                By default, the center of the plume at y=0 is plotted.
            time (float): Point of time for the plot. Will show the closest time step to given value.
                By default, last point in time is plotted.
            legend_names (str | list): List of legend names as strings, in the same order as given models.
                By default, no legend is shown.
            animate (bool): If True, animation of contaminant plume until given time is shown. If multiple models are
                given as input, dt should be the same for each one to ensure accurate animation. Default is False.
            **kwargs : Arguments to be passed to plt.plot().

        """
        pline.centerline(self, y_position=y_position, time=time, legend_names=legend_names, animate=animate, **kwargs)

    def transverse(self, x_position, time=None, legend_names=None, animate=False, **kwargs):
        """Plot concentration distribution as a line horizontal transverse to the plume extent.

        Args:
            model : Model object from mibitrans.transport, or list of model objects.
            x_position : x-position along the plume (longitudinal direction) for the plot.
            time (float): Point of time for the plot. Will show the closest time step to given value.
                By default, last point in time is plotted.
            legend_names (str | list): List of legend names as strings, in the same order as given models.
                By default, no legend is shown.
            animate (bool): If True, animation of contaminant plume until given time is shown. If multiple models are
                given as input, dt should be the same for each one to ensure accurate animation. Default is False.
            **kwargs : Arguments to be passed to plt.plot().
        """
        pline.transverse(self, x_position=x_position, time=time, legend_names=legend_names, animate=animate, **kwargs)

    def breakthrough(self, x_position, y_position=0, legend_names=None, animate=False, **kwargs):
        """Plot contaminant breakthrough curve at given x and y position in model domain.

        Args:
            model : Model object from mibitrans.transport, or list of model objects.
            x_position : x-position along the plume (longitudinal direction).
            y_position : y-position across the plume (transverse horizontal direction).
                By default, at the center of the plume (at y=0).
            legend_names (str | list): List of legend names as strings, in the same order as given models.
                By default, no legend is shown.
            animate (bool): If True, animation of contaminant plume until given time is shown. If multiple models are
                given as input, dt should be the same for each one to ensure accurate animation. Default is False.
            **kwargs : Arguments to be passed to plt.plot().
        """
        pline.breakthrough(
            self, x_position=x_position, y_position=y_position, legend_names=legend_names, animate=animate, **kwargs
        )

    def plume_2d(self, time=None, animate=False, **kwargs):
        """Plot contaminant plume as a 2D colormesh, at a specified time.

        Args:
            model : Model object from mibitrans.transport.
            time (float): Point of time for the plot. Will show the closest time step to given value.
                By default, last point in time is plotted.
            animate (bool): If True, animation of contaminant plume until given time is shown.
            **kwargs : Arguments to be passed to plt.pcolormesh().

        Returns a matrix plot of the input plume as object.
        """
        psurf.plume_2d(self, time=time, animate=animate, **kwargs)

    def plume_3d(self, time=None, animate=False, **kwargs):
        """Plot contaminant plume as a 3D surface, at a specified time.

        Args:
            model : Model object from mibitrans.transport.
            time (float): Point of time for the plot. Will show the closest time step to given value.
                By default, last point in time is plotted.
            animate (bool): If True, animation of contaminant plume until given time is shown.
            **kwargs : Arguments to be passed to plt.plot_surface().

        Returns:
            ax (matplotlib.axes._axes.Axes) : Returns matplotlib axes object of plume plot.
        """
        psurf.plume_3d(self, time=time, animate=animate, **kwargs)


class Domenico(Transport3D, ABC):
    """Parent class that for all analytical solutions using the Domenico (1987) analytical model.

    Domenico, P. A. (1987). An analytical model for multidimensional transport of a decaying contaminant species.
    Journal of Hydrology, 91(1-2), 49-58.
    """

    def __init__(
        self,
        hydrological_parameters,
        attenuation_parameters,
        source_parameters,
        model_parameters,
        verbose=False,
    ):
        """Initialize object and run model.

        Args:
            hydrological_parameters (mibitrans.data.parameters.HydrologicalParameters) : Dataclass object containing
                hydrological parameters from HydrologicalParameters.
            attenuation_parameters (mibitrans.data.read.AttenuationParameters) : Dataclass object containing adsorption,
                degradation and diffusion parameters from AttenuationParameters.
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
        super().__init__(hydrological_parameters, attenuation_parameters, source_parameters, model_parameters, verbose)
        if self._att_pars.diffusion != 0:
            warnings.warn("Domenico model does not consider molecular diffusion.", UserWarning)

    @abstractmethod
    def _calculate_cxyt(self, xxx, yyy, ttt):
        pass

    def run(self):
        """Calculate the concentration for all discretized x, y and t using the analytical transport model."""
        self._pre_run_initialization_parameters()
        self.cxyt = self._calculate_cxyt(self.xxx, self.yyy, self.ttt)

    def _eq_x_term(self, xxx, ttt, decay_sqrt=1):
        return erfc(
            (xxx - self._hyd_pars.velocity * ttt * decay_sqrt)
            / (2 * np.sqrt(self._hyd_pars.alpha_x * self._hyd_pars.velocity * ttt))
        )

    def _eq_additional_x(self, xxx, ttt):
        return np.exp(xxx * self._hyd_pars.velocity / (self._hyd_pars.alpha_x * self._hyd_pars.velocity)) * (
            erfc(
                xxx
                + self._hyd_pars.velocity * ttt / (2 * np.sqrt(self._hyd_pars.alpha_x * self._hyd_pars.velocity * ttt))
            )
        )

    def _eq_z_term(self, xxx):
        inner_term = self._src_pars.depth / (2 * np.sqrt(self._hyd_pars.alpha_z * xxx))
        return erf(inner_term) - erf(-inner_term)

    def _eq_source_decay(self, xxx, ttt):
        term = np.exp(-self.k_source * (ttt - xxx / self._hyd_pars.velocity))
        # Term can be max 1; can not have 'generation' of solute ahead of advection.
        return np.where(term > 1, 1, term)

    def _eq_y_term(self, i, xxx, yyy):
        div_term = 2 * np.sqrt(self._hyd_pars.alpha_y * xxx)
        term = erf((yyy + self.y_source[i]) / div_term) - erf((yyy - self.y_source[i]) / div_term)
        term[np.isnan(term)] = 0
        return term


class Karanovic(Transport3D):
    """Parent class that for all models using the exact analytical solution described in Karanovic (2007).

    Karanovic, M., Neville, C. J., & Andrews, C. B. (2007). BIOSCREEN‐AT: BIOSCREEN with an exact analytical solution.
    Groundwater, 45(2), 242-245.
    """

    def __init__(
        self,
        hydrological_parameters,
        attenuation_parameters,
        source_parameters,
        model_parameters,
        verbose=False,
    ):
        """Initialize object and run model.

        Args:
            hydrological_parameters (mibitrans.data.parameters.HydrologicalParameters) : Dataclass object containing
                hydrological parameters from HydrologicalParameters.
            attenuation_parameters (mibitrans.data.read.AttenuationParameters) : Dataclass object containing adsorption,
                degradation and diffusion parameters from AttenuationParameters.
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
        super().__init__(hydrological_parameters, attenuation_parameters, source_parameters, model_parameters, verbose)

    @abstractmethod
    def _calculate_cxyt(self):
        pass

    def _pre_run_initialization_parameters(self):
        super()._pre_run_initialization_parameters()
        self.disp_x = self._hyd_pars.alpha_x * self.rv + self._att_pars.diffusion
        self.disp_y = self._hyd_pars.alpha_y * self.rv + self._att_pars.diffusion
        self.disp_z = self._hyd_pars.alpha_z * self.rv + self._att_pars.diffusion
        # self.integral_term = np.zeros(self.ttt.shape)
        # Stores integral error for each time step and source zone
        self.error_size = np.zeros((len(self._src_pars.source_zone_boundary), len(self.t)))

    def run(self):
        """Calculate the concentration for all discretized x, y and t using the analytical transport model."""
        self._pre_run_initialization_parameters()
        self.cxyt = self._calculate_cxyt()

    def _eq_integrand(self, t, sz):
        term = 1 / (t ** (3 / 2)) * self._eq_x_exp_term(t) * self._eq_y_term(t, sz) * self._eq_z_term(t)
        term[np.isnan(term)] = 0
        return term

    def _eq_x_exp_term(self, t):
        term = np.exp(
            (-self.k_source - self._att_pars.decay_rate) * t
            - (self.xxx[0, :, 1:] - self.rv * t) ** 2 / (4 * self.disp_x * t)
        )
        term[np.isnan(term)] = 0
        return term

    def _eq_y_term(self, t, sz):
        div_term = 2 * np.sqrt(self.disp_y * t)
        term = erfc((self.yyy[0, :, 1:] - self.y_source[sz]) / div_term) - erfc(
            (self.yyy[0, :, 1:] + self.y_source[sz]) / div_term
        )
        term[np.isnan(term)] = 0
        return term

    def _eq_z_term(self, t):
        if t == 0 or self.disp_z == 0:
            inner_term = 2
        else:
            inner_term = self._src_pars.depth / (2 * np.sqrt(self.disp_z * t))
        return erfc(-inner_term) - erfc(inner_term)

    def _eq_source_term(self, sz):
        return (
            self.c_source[sz]
            * self.xxx[:, :, 1:]
            / (8 * np.sqrt(np.pi * self.disp_x))
            * np.exp(-self.k_source * self.ttt[:, :, 1:])
        )

    def _eq_source_zero(self, sz):
        zone_location = np.where(abs(self.yyy[:, :, 0]) <= self.y_source[sz], 1, 0)
        return self.c_source[sz] * zone_location * np.exp(-self.k_source * self.ttt[:, :, 0])

    def _eq_integral_term(self, sz):
        integral_term = np.zeros(self.ttt.shape)
        for j in range(len(self.t)):
            if self.verbose:
                print("integrating for t =", self.t[j], "days")
            if j == 0:
                lower_bound = 0
            else:
                lower_bound = self.t[j - 1]
            upper_bound = self.t[j]
            integral_term[j, :, 1:], self.error_size[sz, j] = quad_vec(
                self._eq_integrand, lower_bound, upper_bound, limit=10000 / len(self.t), args=(sz,)
            )
        integral_sum = np.cumsum(integral_term, axis=0)
        return integral_sum

    def _eq_superposition(self):
        cxyt = self.cxyt.copy()
        for sz in range(len(self.c_source)):
            if self.verbose:
                print("integrating for source zone ", sz)
            integral_sum = self._eq_integral_term(sz)
            source_term = self._eq_source_term(sz)
            cxyt[:, :, 1:] += integral_sum[:, :, 1:] * source_term
            # If x=0, equation resolves to c=0, therefore, x=0 needs to be evaluated separately
            cxyt[:, :, 0] += self._eq_source_zero(sz)
        return cxyt

    def sample(self, x_position, y_position, time):
        """Give concentration at any given position and point in time.

        Args:
            x_position (float): x position in domain extent [m].
            y_position (float): y position in domain extent [m].
            time (float): time for which concentration is sampled [days].

        Returns:
            concentration (float): concentration at given position and point in time [g/m^3].

        """
        # Different sample method than parent class, as calculations from _calculate use array indices
        for par, value in locals().items():
            if par != "self":
                validate_input_values(par, value)

        if not self.has_run and not self.initialized:
            self._pre_run_initialization_parameters()

        def integrand(t, sz):
            div_term = 2 * np.sqrt(self.disp_y * t**4)
            inner_term = self._src_pars.depth / (2 * np.sqrt(self.disp_z * t**4))
            integrand_results = (
                1
                / (t**3)
                * (
                    np.exp(
                        (-self.k_source - self._att_pars.decay_rate) * t**4
                        - (x_position - self.rv * t**4) ** 2 / (4 * self.disp_x * t**4)
                    )
                    * (
                        erfc((y_position - self.y_source[sz]) / div_term)
                        - erfc((y_position + self.y_source[sz]) / div_term)
                    )
                    * (erfc(-inner_term) - erfc(inner_term))
                )
            )
            return integrand_results

        conc_array = np.zeros(len(self.c_source))
        error_array = np.zeros(len(self.c_source))
        time = time ** (1 / 4)
        with np.errstate(divide="ignore", invalid="ignore"):
            for sz in range(len(self.c_source)):
                integral_term, error = quad(integrand, 0, time, limit=10000, args=(sz,))
                source_term = (
                    self.c_source[sz] * x_position / (8 * np.sqrt(np.pi * self.disp_x)) * np.exp(-self.k_source * time)
                )
                conc_array[sz] = 4 * integral_term * source_term
                error_array[sz] = error
            concentration = np.sum(conc_array)
        return concentration
