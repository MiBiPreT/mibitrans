import copy
import warnings
from abc import ABC
from abc import abstractmethod
import numpy as np
import mibitrans.data.parameters
from mibitrans.data.parameter_information import ElectronAcceptors
from mibitrans.data.parameter_information import UtilizationFactor
from mibitrans.data.parameter_information import util_to_conc_name
from mibitrans.visualize import plot_line as pline
from mibitrans.visualize import plot_surface as psurf


class Transport3D(ABC):
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
        self._decay_rate = self._att_pars.decay_rate

        self.verbose = verbose

        self._observe_input_dataclass_change()

        self.has_run = False
        self.initialized = False
        self._mode = "linear"
        self._electron_acceptors = None
        self._utilization_factor = None
        self.biodegradation_capacity = None
        self.cxyt_noBC = None
        self._pre_run_initialization_parameters()

    @property
    def input_parameters(self):
        """Return the input arguments for the model in the form of a dictionary, based on current values."""
        return dict(
            hydrological_parameters=self.hydrological_parameters,
            attenuation_parameters=self.attenuation_parameters,
            source_parameters=self.source_parameters,
            model_parameters=self.model_parameters,
            verbose=self.verbose,
        )

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

    @property
    def mode(self):
        """Model mode property. Either 'linear' or 'instant_reaction'."""
        return self._mode

    @mode.setter
    def mode(self, value):
        match value:
            case "linear" | "linear decay" | "linear_decay" | 0:
                self._mode = "linear"
                self.initialized = False
            case "instant" | "instant_reaction" | "instant reaction" | 1:
                if self._electron_acceptors is None or self._utilization_factor is None:
                    raise ValueError(
                        "Model mode was set to 'instant reaction', but electron acceptor parameters are "
                        "missing. Use the instant_reaction method to supply the electron acceptor "
                        "concentrations."
                    )
                self._mode = "instant_reaction"
                self.initialized = False
            case _:
                warnings.warn(f"Mode '{value}' not recognized. Defaulting to 'linear' instead.", UserWarning)
                self._mode = "linear"
                self.initialized = False

    @property
    def electron_acceptors(self):
        """Return dictionary of electron acceptor parameters."""
        return self._electron_acceptors.dictionary

    @property
    def utilization_factor(self):
        """Return dictionary of utilization factor property."""
        return self._utilization_factor.dictionary

    @property
    def relative_cxyt(self):
        """Compute relative concentration c(x,y,t)/c0, where c0 is the maximum source zone concentration at t=0."""
        maximum_concentration = np.max(self.source_parameters.source_zone_concentration)
        relative_cxyt = self.cxyt / maximum_concentration
        return relative_cxyt

    @property
    @abstractmethod
    def short_description(self):
        """Short string describing model type."""
        pass

    @abstractmethod
    def run(self):
        """Method that runs the model and ensures that initialisation is performed."""
        pass

    @abstractmethod
    def sample(self, x_position, y_position, t_position):
        """Method that calculates concentration at single, specified location in model domain."""
        pass

    @abstractmethod
    def _calculate_concentration_for_all_xyt(self) -> np.ndarray:
        """Method that calculates and return concentration array for all model x, y and t."""
        pass

    def _observe_input_dataclass_change(self):
        """Keeps track of input dataclass changes, and ensures re-initialization for following model runs."""
        self.hydrological_parameters._on_change = lambda: self._check_and_reset_when_input_dataclass_change(
            "hydrological_parameters", self._hyd_pars
        )
        self.attenuation_parameters._on_change = lambda: self._check_and_reset_when_input_dataclass_change(
            "attenuation_parameters", self._att_pars
        )
        self.source_parameters._on_change = lambda: self._check_and_reset_when_input_dataclass_change(
            "source_parameters", self._src_pars
        )
        self.model_parameters._on_change = lambda: self._check_and_reset_when_input_dataclass_change(
            "model_parameters", self._mod_pars
        )

    def _check_and_reset_when_input_dataclass_change(self, key, value):
        """Remove output and unflag initialization when input dataclass changes are observed."""
        self._check_input_dataclasses(key, value)
        if self.verbose:
            print(key, " got changed")
        self.initialized = False
        if self.has_run:
            self.cxyt = np.zeros((len(self.t), len(self.y), len(self.x)))
            if self.cxyt_noBC is not None:
                self.cxyt_noBC = 0
            self.has_run = False
            if self.verbose:
                print(f"Parameter '{key}' has changed — resetting output")

    def _pre_run_initialization_parameters(self):
        """Parameter initialization for model."""
        # One-dimensional model domain arrays
        self.x = np.arange(0, self._mod_pars.model_length + self._mod_pars.dx, self._mod_pars.dx)
        self.y = self._calculate_y_discretization()
        self.t = np.arange(self._mod_pars.dt, self._mod_pars.model_time + self._mod_pars.dt, self._mod_pars.dt)

        # Three-dimensional model domain arrays
        self.xxx = self.x[None, None, :]
        self.yyy = self.y[None, :, None]
        self.ttt = self.t[:, None, None]

        self.rv = self._hyd_pars.velocity / self._att_pars.retardation

        # cxyt is concentration output array
        self.cxyt = np.zeros((len(self.t), len(self.y), len(self.x)))

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
        if self._mode == "instant_reaction":
            self.c_source[-1] += self.biodegradation_capacity
            self._decay_rate = 0
        else:
            self._decay_rate = self._att_pars.decay_rate

        self.initialized = True

    def _calculate_source_decay(self):
        """Calculate source decay/depletion."""
        if isinstance(self._src_pars.total_mass, (float, int)):
            Q, c0_avg = self._calculate_discharge_and_average_source_zone_concentration()
            k_source = Q * c0_avg / self._src_pars.total_mass
        # If source mass is not a float, it is an infinite source, therefore, no source decay takes place.
        else:
            k_source = 0

        return k_source

    def _calculate_discharge_and_average_source_zone_concentration(self):
        """Calculate the average source zone concentration, and discharge through source zone."""
        if self._mode == "instant_reaction":
            bc = self.biodegradation_capacity
        else:
            bc = 0
        y_src = np.zeros(len(self._src_pars.source_zone_boundary) + 1)
        y_src[1:] = self._src_pars.source_zone_boundary
        c_src = self._src_pars.source_zone_concentration
        Q = self._hyd_pars.velocity * self._hyd_pars.porosity * self._src_pars.depth * np.max(y_src) * 2

        weighted_conc = np.zeros(len(self._src_pars.source_zone_boundary))
        for i in range(len(self._src_pars.source_zone_boundary)):
            weighted_conc[i] = (y_src[i + 1] - y_src[i]) * c_src[i]

        c0_avg = bc + np.sum(weighted_conc) / np.max(y_src)

        return Q, c0_avg

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

    def _calculate_biodegradation_capacity(self):
        """Determine biodegradation capacity based on electron acceptor concentrations and utilization factor."""
        biodegradation_capacity = 0
        for key, item in self._utilization_factor.dictionary.items():
            biodegradation_capacity += getattr(self._electron_acceptors, util_to_conc_name[key]) / item

        return biodegradation_capacity

    def instant_reaction(
        self,
        electron_acceptors: list | np.ndarray | dict | ElectronAcceptors,
        utilization_factor: list | np.ndarray | dict | UtilizationFactor = UtilizationFactor(
            util_oxygen=3.14, util_nitrate=4.9, util_ferrous_iron=21.8, util_sulfate=4.7, util_methane=0.78
        ),
    ):
        """Enable and set up parameters for instant reaction model.

        Instant reaction model assumes that biodegradation is an instantaneous process compared to the groundwater flow
        velocity. The biodegradation is assumed to be governed by the availability of electron acceptors, and quantified
        using  stoichiometric relations from the degradation reactions. Considered are concentrations of acceptors
        Oxygen, Nitrate and Sulfate, and reduced species Ferrous Iron and Methane.

        Args:
            electron_acceptors (ElectronAcceptors): ElectronAcceptor dataclass containing electron acceptor
                concentrations. Alternatively provided as list, numpy array or dictionary corresponding with
                delta_oxygen, delta_nitrate, ferrous_iron, delta_sulfate and methane. For more information, see
                documentation for ElectronAcceptors.
            utilization_factor (UtilizationFactor, optional): UtilizationFactor dataclass containing electron acceptor
                utilization factors. Alternatively provided as list, numpy array or dictionary corresponding with
                information, see documentation of UtilizationFactor. By default, electron acceptor utilization factors
                for a BTEX mixture are used, based on values by Wiedemeier et al. (1995).
        """
        self._electron_acceptors, self._utilization_factor = _check_instant_reaction_acceptor_input(
            electron_acceptors, utilization_factor
        )
        self._mode = "instant_reaction"
        self.biodegradation_capacity = self._calculate_biodegradation_capacity()
        self.cxyt_noBC = 0
        self._pre_run_initialization_parameters()

    def _check_model_mode_before_run(self):
        # Reset concentration array to make sure it is empty before calculation.
        self.cxyt = np.zeros((len(self.t), len(self.y), len(self.x)))
        if not self.initialized:
            self._pre_run_initialization_parameters()
        if self._mode == "linear":
            if self.biodegradation_capacity is not None:
                warnings.warn(
                    "Instant reaction parameters are present while model mode is linear. "
                    "Make sure that this is indeed the desired model."
                )
        if self._mode == "instant_reaction":
            if self.biodegradation_capacity is None:
                raise ValueError(
                    "Instant reaction parameters are not present. "
                    "Please provide them with the 'instant_reaction' class method."
                )

    def centerline(self, y_position=0, time=None, relative_concentration=False, animate=False, **kwargs):
        """Plot center of contaminant plume of this model, at a specified time and y position.

        Args:
            y_position (float, optional): y-position across the plume (transverse horizontal direction) for the plot.
                By default, the center of the plume at y=0 is plotted.
            time (float, optional): Point of time for the plot. Will show the closest time step to given value.
                By default, last point in time is plotted.
            relative_concentration (bool, optional) : If set to True, will plot concentrations relative to maximum
                source zone concentrations at t=0. By default, absolute concentrations are shown.
            animate (bool, optional): If True, animation of contaminant plume until given time is shown. Default is
                False.
            **kwargs : Arguments to be passed to plt.plot().

        """
        if animate:
            anim = pline.centerline(
                self,
                y_position=y_position,
                time=time,
                relative_concentration=relative_concentration,
                animate=animate,
                **kwargs,
            )
            return anim
        else:
            pline.centerline(
                self,
                y_position=y_position,
                time=time,
                relative_concentration=relative_concentration,
                animate=animate,
                **kwargs,
            )
            return None

    def transverse(self, x_position, time=None, relative_concentration=False, animate=False, **kwargs):
        """Plot concentration distribution as a line horizontal transverse to the plume extent.

        Args:
            x_position : x-position along the plume (longitudinal direction) for the plot.
            time (float): Point of time for the plot. Will show the closest time step to given value.
                By default, last point in time is plotted.
            relative_concentration (bool, optional) : If set to True, will plot concentrations relative to maximum
                source zone concentrations at t=0. By default, absolute concentrations are shown.
            animate (bool, optional): If True, animation of contaminant plume until given time is shown. Default is
                False.
            **kwargs : Arguments to be passed to plt.plot().
        """
        if animate:
            anim = pline.transverse(
                self,
                x_position=x_position,
                time=time,
                relative_concentration=relative_concentration,
                animate=animate,
                **kwargs,
            )
            return anim
        else:
            pline.transverse(
                self,
                x_position=x_position,
                time=time,
                relative_concentration=relative_concentration,
                animate=animate,
                **kwargs,
            )
            return None

    def breakthrough(self, x_position, y_position=0, relative_concentration=False, animate=False, **kwargs):
        """Plot contaminant breakthrough curve at given x and y position in model domain.

        Args:
            x_position : x-position along the plume (longitudinal direction).
            y_position : y-position across the plume (transverse horizontal direction).
                By default, at the center of the plume (at y=0).
            relative_concentration (bool, optional) : If set to True, will plot concentrations relative to maximum
                source zone concentrations at t=0. By default, absolute concentrations are shown.
            animate (bool, optional): If True, animation of contaminant plume until given time is shown. Default is
                False.
            **kwargs : Arguments to be passed to plt.plot().
        """
        if animate:
            anim = pline.breakthrough(
                self,
                x_position=x_position,
                y_position=y_position,
                relative_concentration=relative_concentration,
                animate=animate,
                **kwargs,
            )
            return anim
        else:
            pline.breakthrough(
                self,
                x_position=x_position,
                y_position=y_position,
                relative_concentration=relative_concentration,
                animate=animate,
                **kwargs,
            )
            return None

    def plume_2d(self, time=None, relative_concentration=False, animate=False, **kwargs):
        """Plot contaminant plume as a 2D colormesh, at a specified time.

        Args:
            time (float): Point of time for the plot. Will show the closest time step to given value.
                By default, last point in time is plotted.
            relative_concentration (bool, optional) : If set to True, will plot concentrations relative to maximum
                source zone concentrations at t=0. By default, absolute concentrations are shown.
            animate (bool, optional): If True, animation of contaminant plume until given time is shown. Default is
                False.
            **kwargs : Arguments to be passed to plt.pcolormesh().

        Returns a matrix plot of the input plume as object.
        """
        anim = psurf.plume_2d(self, time=time, relative_concentration=relative_concentration, animate=animate, **kwargs)
        return anim

    def plume_3d(self, time=None, relative_concentration=False, animate=False, **kwargs):
        """Plot contaminant plume as a 3D surface, at a specified time.

        Args:
            time (float): Point of time for the plot. Will show the closest time step to given value.
                By default, last point in time is plotted.
            relative_concentration (bool, optional) : If set to True, will plot concentrations relative to maximum
                source zone concentrations at t=0. By default, absolute concentrations are shown.
            animate (bool, optional): If True, animation of contaminant plume until given time is shown. Default is
                False.
            **kwargs : Arguments to be passed to plt.plot_surface().

        Returns:
            ax (matplotlib.axes._axes.Axes) : Matplotlib Axes object of plume plot.
                or if animate == True
            anim (matplotib.animation.FuncAnimation) : Matplotlib FuncAnimation object of plume plot.
        """
        ax_or_anim = psurf.plume_3d(
            self, time=time, relative_concentration=relative_concentration, animate=animate, **kwargs
        )
        return ax_or_anim


def _check_instant_reaction_acceptor_input(electron_acceptors, utilization_factor):
    if isinstance(electron_acceptors, (list, np.ndarray)):
        electron_acceptors_out = ElectronAcceptors(*electron_acceptors)
    elif isinstance(electron_acceptors, dict):
        electron_acceptors_out = ElectronAcceptors(**electron_acceptors)
    elif isinstance(electron_acceptors, mibitrans.data.parameter_information.ElectronAcceptors):
        electron_acceptors_out = electron_acceptors
    else:
        raise TypeError(
            f"electron_acceptors must be a list, dictionary or ElectronAcceptors dataclass, but is "
            f"{type(electron_acceptors)} instead."
        )

    if isinstance(utilization_factor, (list, np.ndarray)):
        utilization_factor_out = UtilizationFactor(*utilization_factor)
    elif isinstance(utilization_factor, dict):
        utilization_factor_out = UtilizationFactor(**utilization_factor)
    elif isinstance(utilization_factor, mibitrans.data.parameter_information.UtilizationFactor):
        utilization_factor_out = utilization_factor
    else:
        raise TypeError(
            f"utilization_factor must be a list, dictionary or UtilizationFactor dataclass, but is "
            f"{type(utilization_factor)} instead."
        )

    return electron_acceptors_out, utilization_factor_out
