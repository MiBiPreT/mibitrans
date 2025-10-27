import numpy as np
from mibitrans.transport.model_parent import Karanovic


class LinearDecay(Karanovic):
    """Calculate contaminant transport with linear decay using the exact analytical solution."""

    def __init__(
        self,
        hydrological_parameters,
        attenuation_parameters,
        source_parameters,
        model_parameters,
        auto_run=True,
        verbose=False,
    ):
        """Initialize object and run model.

        Args:
            hydrological_parameters (mibitrans.data.parameters.HydrologicalParameters) : Dataclass object containing
                hydrological parameters from HydrologicalParameters.
            attenuation_parameters (mibitrans.data.parameters.AttenuationParameters) : Dataclass object containing
                adsorption, degradation and diffusion parameters from AttenuationParameters.
            source_parameters (mibitrans.data.parameters.SourceParameters) : Dataclass object containing source
                parameters from SourceParameters.
            model_parameters (mibitrans.data.parameters.ModelParameters) : Dataclass object containing model parameters
                from ModelParameters.
            auto_run (bool, optional): Automatically run model when initialized. Defaults to True.
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
            sample : Give concentration at any given position and point in time.

        Raises:
            TypeError : If input is not of the correct Dataclass.

        """
        super().__init__(hydrological_parameters, attenuation_parameters, source_parameters, model_parameters, verbose)
        if auto_run:
            self.cxyt = self._calculate_concentration_for_all_xyt()

    @property
    def short_description(self):
        """Short string describing model type."""
        return "Karanovic Linear Decay"

    def _calculate_concentration_for_all_xyt(self):
        with np.errstate(divide="ignore", invalid="ignore"):
            cxyt = self._equation_source_superposition()
        self.has_run = True
        return cxyt


class NoDecay(Karanovic):
    """Calculate contaminant transport with no decay using the exact analytical solution."""

    def __init__(
        self,
        hydrological_parameters,
        attenuation_parameters,
        source_parameters,
        model_parameters,
        auto_run=True,
        verbose=False,
    ):
        """Initialize object and run model.

        Args:
            hydrological_parameters (mibitrans.data.parameters.HydrologicalParameters) : Dataclass object containing
            hydrological parameters from HydrologicalParameters.
            attenuation_parameters (mibitrans.data.parameters.AttenuationParameters) : Dataclass object containing
                adsorption, degradation and diffusion parameters from AttenuationParameters.
            source_parameters (mibitrans.data.parameters.SourceParameters) : Dataclass object containing source
                parameters from SourceParameters.
            model_parameters (mibitrans.data.parameters.ModelParameters) : Dataclass object containing model parameters
                from ModelParameters.
            auto_run (bool, optional): Automatically run model when initialized. Defaults to True.
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
            sample : Give concentration at any given position and point in time.

        Raises:
            TypeError : If input is not of the correct Dataclass.

        """
        # attenuation_parameters = copy.copy(attenuation_parameters)
        # attenuation_parameters.decay_rate = 0
        # super().__init__(
        #     hydrological_parameters, attenuation_parameters, source_parameters, model_parameters, auto_run, verbose
        # )
        super().__init__(hydrological_parameters, attenuation_parameters, source_parameters, model_parameters, verbose)
        self._att_pars.decay_rate = 0
        if auto_run:
            self.cxyt = self._calculate_concentration_for_all_xyt()

    @property
    def short_description(self):
        """Short string describing model type."""
        return "Karanovic No Decay"

    def _calculate_concentration_for_all_xyt(self):
        with np.errstate(divide="ignore", invalid="ignore"):
            cxyt = self._equation_source_superposition()
        self.has_run = True
        return cxyt


class InstantReaction(Karanovic):
    """Calculate contaminant transport with instant reaction degradation using the exact analytical solution."""

    def __init__(
        self,
        hydrological_parameters,
        attenuation_parameters,
        source_parameters,
        model_parameters,
        auto_run=True,
        verbose=False,
    ):
        """Initialize object and run model.

        Args:
            hydrological_parameters (mibitrans.data.parameters.HydrologicalParameters) : Dataclass object containing
            hydrological parameters from HydrologicalParameters.
            attenuation_parameters (mibitrans.data.parameters.AttenuationParameters) : Dataclass object containing
                adsorption, degradation and diffusion parameters from AttenuationParameters.
            source_parameters (mibitrans.data.parameters.SourceParameters) : Dataclass object containing source
                parameters from SourceParameters.
            model_parameters (mibitrans.data.parameters.ModelParameters) : Dataclass object containing model parameters
                from ModelParameters.
            auto_run (bool, optional): Automatically run model when initialized. Defaults to True.
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
        if auto_run:
            self.cxyt = self._calculate_concentration_for_all_xyt()

    @property
    def short_description(self):
        """Short string describing model type."""
        return "Karanovic Instant Reaction"

    def _instant_initialization(self):
        self._att_pars._require_electron_acceptor()
        self._att_pars.decay_rate = 0
        self.biodegradation_capacity = self._calculate_biodegradation_capacity()
        self.k_source = self._calculate_source_decay(self.biodegradation_capacity)
        # Source decay calculation uses self.c_source, therefore, addition of biodegradation_capacity to
        # outer source zone after calculation of k_source
        self.c_source[-1] += self.biodegradation_capacity
        self.cxyt_noBC = 0

    def _pre_run_initialization_parameters(self):
        super()._pre_run_initialization_parameters()
        self._instant_initialization()

    def _calculate_concentration_for_all_xyt(self):
        with np.errstate(divide="ignore", invalid="ignore"):
            cxyt = self._equation_source_superposition()
            self.cxyt_noBC = cxyt.copy()
            cxyt -= self.biodegradation_capacity
            cxyt = np.where(cxyt < 0, 0, cxyt)
        self.has_run = True
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
        concentration_noBC = super().sample(x_position, y_position, time)
        concentration = concentration_noBC - self.biodegradation_capacity
        if concentration < 0:
            concentration = 0
        return concentration
