import copy
import warnings
from abc import ABC
from abc import abstractmethod
import numpy as np
import mibitrans
import mibitrans.data.parameters
from mibitrans.analysis.mass_balance import MassBalance
from mibitrans.analysis.parameter_calculations import calculate_discharge_and_average_source_zone_concentration
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
        self._check_input_dataclasses("hydrological_parameters", value)
        self._hyd_pars = copy.copy(value)

    @property
    def attenuation_parameters(self):
        """Rename to shorthand form of attenuation_parameters inside class for ease of use."""
        return self._att_pars

    @attenuation_parameters.setter
    def attenuation_parameters(self, value):
        self._check_input_dataclasses("attenuation_parameters", value)
        self._att_pars = copy.copy(value)

    @property
    def source_parameters(self):
        """Rename to shorthand form of source_parameters inside class for ease of use."""
        return self._src_pars

    @source_parameters.setter
    def source_parameters(self, value):
        self._check_input_dataclasses("source_parameters", value)
        self._src_pars = copy.copy(value)

    @property
    def model_parameters(self):
        """Rename to shorthand form of model_parameters inside class for ease of use."""
        return self._mod_pars

    @model_parameters.setter
    def model_parameters(self, value):
        self._check_input_dataclasses("model_parameters", value)
        self._mod_pars = copy.copy(value)

    @property
    def mode(self):
        """Model mode property. Either 'linear' or 'instant_reaction'."""
        return self._mode

    @mode.setter
    def mode(self, value):
        match value:
            case "linear" | "linear decay" | "linear_decay" | 0:
                self._mode = "linear"
            case "instant" | "instant_reaction" | "instant reaction" | 1:
                if self._electron_acceptors is None or self._utilization_factor is None:
                    raise ValueError(
                        "Model mode was set to 'instant reaction', but electron acceptor parameters are "
                        "missing. Use the instant_reaction method to supply the electron acceptor "
                        "concentrations."
                    )
                self._mode = "instant_reaction"
            case _:
                warnings.warn(f"Mode '{value}' not recognized. Defaulting to 'linear' instead.", UserWarning)
                self._mode = "linear"

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

        if (
            self._att_pars.bulk_density is not None
            and self._att_pars.partition_coefficient is not None
            and self._att_pars.fraction_organic_carbon is not None
        ):
            self._att_pars.calculate_retardation(self._hyd_pars.porosity)

        self.rv = self._hyd_pars.velocity / self._att_pars.retardation

        # cxyt is concentration output array
        self.cxyt = np.zeros((len(self.t), len(self.y), len(self.x)))

        # Calculate retardation if not already specified in adsorption_parameters
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

    def _calculate_source_decay(self):
        """Calculate source decay/depletion."""
        if self._src_pars.total_mass != np.inf:
            Q, c0_avg = calculate_discharge_and_average_source_zone_concentration(self)
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


class Results:
    """Object that holds model results and input parameters for individual runs."""

    def __init__(self, model):
        """Records input parameters and resulting output based given model.

        Args:
            model (Transport3D): Model object from which to initialize results. Should be child class of Transport3D.

        Properties:
            model_type (Transport3D) : Class instance of model used to generate results.
            short_description (str) : Short description of model.
            x (np.ndarray) : Numpy array with model x (longitudinal direction) discretization, corresponding to
                model_parameters, with distance in [m].
            y (np.ndarray) : Numpy array with model y (transverse horizontal direction) discretization, corresponding to
                model_parameters, with distance in [m].
            t (np.ndarray) : Numpy array with model t (time) discretization, corresponding to
                model_parameters, with time in [days].
            hydrological_parameters (HydrologicalParameters) : Dataclass holding the hydrological parameters used to
                run the model.
            attenuation_parameters (AttenuationParameters) : Dataclass holding the attenuation parameters used to run
                the model.
            source_parameters (SourceParameters) : Dataclass holding the source parameters used to run the model.
            model_parameters (ModelParameters): Dataclass holding the model parameters used to run the model.
            electron_acceptors (ElectronAcceptors): Dataclass holding the electron acceptor concentrations used to run
                the model. Only for instant reaction, None for other models.
            utilization_factor (UtilizationFactor): Dataclass holding the electron acceptor utilization factors used to
                run the model. Only for instant reaction, None for other models.
            mode (str) : Model mode of the used model. Either 'linear' or 'instant_reaction'
            rv (float) : Retarded flow velocity, as v / R [m/day].
            k_source (float) : Source depletion rate [1/days]. For infinite source mass, k_source = 0, and therefore, no
                source depletion takes place.
            c_source (np.ndarray) : Initial nett source zone concentrations. For multiple source zones, nett
                concentration in nth source zone is original concentration minus concentration in source zone n - 1. For
                instant reaction model, the biodegradation capacity is added to the outermost source zone.
            biodegradation_capacity (float) : Maximum capacity of biodegradation taking place, based on electron
                acceptor concentrations and utilization factor.
            cxyt (np.ndarray) : Three-dimensional numpy array with concentrations for all x, y and t positions. Indexed
             as cxyt[t,y,x]. In [g/m3].
            relative_cxyt (np.ndarray) : Three-dimensional numpy array with relative concentrations for all x, y and t
                positions. Compared to maximum source zone concentrations.
            cxyt_noBC (np.ndarray) : Three-dimensional numpy array with concentrations for all x, y and t of instant
                reaction models, without subtracting the biodegradation capacity, in [g/m3].
            input_parameters (dict) : Dictionary of input parameter dataclasses for the model. Does not include instant
                reaction parameters.

        Methods:
            centerline : Plot center of contaminant plume, at a specified time and y position.
            transverse : Plot concentration distribution as a line horizontal transverse to the plume extent.
            breakthrough : Plot contaminant breakthrough curve at given x and y position in model domain.
            plume_2d : Plot contaminant plume as a 2D colormesh, at a specified time.
            plume_3d : Plot contaminant plume as a 3D surface, at a specified time.
            mass_balance : Return a mass balance object with source and plume characteristics at given time(s).
        """
        self._model_type = model.__class__
        self._short_description = model.short_description
        self._x = model.x
        self._y = model.y
        self._t = model.t

        # All properties of Transport3D that are objects should be copied; if not copied, changing them in the class
        # object where they originated will also change them here, which is not the intended behaviour.
        self._hydrological_parameters = copy.copy(model.hydrological_parameters)
        self._attenuation_parameters = copy.copy(model.attenuation_parameters)
        self._source_parameters = copy.copy(model.source_parameters)
        self._model_parameters = copy.copy(model.model_parameters)
        self._electron_acceptors = copy.copy(model._electron_acceptors)
        self._utilization_factor = copy.copy(model._utilization_factor)

        self._mode = model.mode
        self._rv = model.rv
        self._k_source = model.k_source
        self._c_source = model.c_source
        self._biodegradation_capacity = model.biodegradation_capacity

        self._cxyt = model.cxyt
        self._relative_cxyt = model.relative_cxyt
        self._cxyt_noBC = model.cxyt_noBC

    @property
    def model_type(self):
        """Class object of the model that generated the results."""
        return self._model_type

    @property
    def short_description(self):
        """Short description of the model that generated the results."""
        return self._short_description

    @property
    def x(self):
        """Model x discretization array."""
        return self._x

    @property
    def y(self):
        """Model y discretization array."""
        return self._y

    @property
    def t(self):
        """Model t discretization array."""
        return self._t

    @property
    def hydrological_parameters(self):
        """Hydrological parameters of the model used for the results."""
        return self._hydrological_parameters

    @property
    def attenuation_parameters(self):
        """Attenuation parameters of the model used for the results."""
        return self._attenuation_parameters

    @property
    def source_parameters(self):
        """Source parameters of the model used for the results."""
        return self._source_parameters

    @property
    def model_parameters(self):
        """Space-time discretization parameters of the model used for the results."""
        return self._model_parameters

    @property
    def electron_acceptors(self):
        """Electron acceptor/byproduct concentrations of the model used for the results."""
        return self._electron_acceptors

    @property
    def utilization_factor(self):
        """Utilization factor of the model used for the results."""
        return self._utilization_factor

    @property
    def mode(self):
        """Model mode used for running the model."""
        return self._mode

    @property
    def rv(self):
        """Retarded flow velocity used in the model."""
        return self._rv

    @property
    def k_source(self):
        """Source depletion rate used in the model."""
        return self._k_source

    @property
    def c_source(self):
        """Nett source zone concentration used in the model."""
        return self._c_source

    @property
    def biodegradation_capacity(self):
        """Biodegradation capacity of the model used for the results. Only for instant reaction models."""
        return self._biodegradation_capacity

    @property
    def cxyt(self):
        """Modelled concentration for all x, y and t, using the input parameters present in this object."""
        return self._cxyt

    @property
    def relative_cxyt(self):
        """Modelled concentration for all x, y and t, divided by the maximum source zone concentration."""
        return self._relative_cxyt

    @property
    def cxyt_noBC(self):
        """Concentration in domain without subtracting biodegradation capacity, in the instant reaction model."""
        return self._cxyt_noBC

    @property
    def input_parameters(self):
        """Return the input arguments for the model in the form of a dictionary, based on current values."""
        return dict(
            hydrological_parameters=self.hydrological_parameters,
            attenuation_parameters=self.attenuation_parameters,
            source_parameters=self.source_parameters,
            model_parameters=self.model_parameters,
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

    def mass_balance(self, time="all", verbose=False):
        """Return a mass balance object with source and plume characteristics at given time(s).

        Args:
            time (float | str): Time at which to initially calculate the mass balance. Either as a value between 0 and
                model end time. Or as 'all', which will calculate mass balance attributes for each time step as arrays.
            verbose (bool, optional): Verbose mode. Defaults to False.

        Returns:
            mass_balance_object: Object of the MassBalance class. Output is accessed through object properties. Can be
                called to change the time of the mass balance.

        The mass balance object has the following properties:
            plume_mass: Mass of the contaminant plume inside the model extent, at the given time(s), in [g].
            source_mass: Mass of the contaminant source at the given time(s), in [g]. No values are given for models
                with infinite source mass.
            delta_source: Difference in mass between contaminant source at given time and source at t = 0, in [g].
            degraded_mass: Mass of plume contaminant degradation at the given time(s), compared to a model without
                degradation, in [g]. Has no value if model does not consider degradation.
            model_without_degradation: Object of model without degradation. Has no value if model does not consider
                degradation.
            instant_reaction_degraded_mass(self): Difference in plume mass instant reaction with and without
                biodegradation capacity subtracted, in [g].
            electron_acceptor_change(self): Change in electron acceptor/byproduct masses at the given time(s), in [g].
                Only for instant reaction.
        """
        return MassBalance(self, time, verbose)


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
