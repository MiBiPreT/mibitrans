"""Author: Jorrit Bakker.

Module calculating the solution to the Domenico (1987) analytical model adapted in BIOSCREEN, for different scenarios.

Domenico, P. A. (1987). An analytical model for multidimensional transport of a decaying contaminant species.
Journal of Hydrology, 91(1-2), 49-58.
"""

import numpy as np
from mibitrans.data.parameter_information import util_to_conc_name
from mibitrans.transport.model_parent import Domenico


class NoDecay(Domenico):
    """Calculate contaminant transport using the Domenico (1987) analytical model without degradation."""

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
            self.cxyt = self._calculate_concentration_for_all_xyt(self.xxx, self.yyy, self.ttt)

    def _calculate_concentration_for_all_xyt(self, xxx, yyy, ttt):
        cxyt = 0
        with np.errstate(divide="ignore", invalid="ignore"):
            x_term = self._equation_term_x(xxx, ttt)
            additional_x = self._equation_term_additional_x(xxx, ttt)
            z_term = self._equation_term_z(xxx)
            source_decay = self._equation_term_source_decay(xxx, ttt)
            for i in range(len(self.c_source)):
                y_term = self._equation_term_y(i, xxx, yyy)
                cxyt_step = 1 / 8 * self.c_source[i] * source_decay * (x_term + additional_x) * y_term * z_term
                cxyt += cxyt_step
        self.has_run = True
        return cxyt


class LinearDecay(Domenico):
    """Calculate contaminant transport using the Domenico (1987) analytical model with a linear decay isotherm."""

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
        self._att_pars._require_linear_decay()
        if auto_run:
            self.cxyt = self._calculate_concentration_for_all_xyt(self.xxx, self.yyy, self.ttt)

    def _calculate_concentration_for_all_xyt(self, xxx, yyy, ttt):
        cxyt = 0
        with np.errstate(divide="ignore", invalid="ignore"):
            decay_sqrt = np.sqrt(1 + 4 * self._att_pars.decay_rate * self._hyd_pars.alpha_x / self._hyd_pars.velocity)
            decay_term = np.exp(xxx * (1 - decay_sqrt) / (self._hyd_pars.alpha_x * 2))
            x_term = self._equation_term_x(xxx, ttt, decay_sqrt)
            z_term = self._equation_term_z(xxx)
            source_decay = self._equation_term_source_decay(xxx, ttt)
            for i in range(len(self.c_source)):
                y_term = self._equation_term_y(i, xxx, yyy)
                cxyt_step = 1 / 8 * self.c_source[i] * source_decay * decay_term * x_term * y_term * z_term
                cxyt += cxyt_step
        self.has_run = True
        return cxyt


class InstantReaction(Domenico):
    """Calculate contaminant transport using the Domenico (1987) analytical model instant reaction biodegradation."""

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
        hydrological_parameters (mibitrans.data.read.HydrologicalParameters) : Dataclass object containing hydrological
            parameters from HydrologicalParameters.
        attenuation_parameters (mibitrans.data.read.AttenuationParameters) : Dataclass object containing adsorption,
                degradation and diffusion parameters from AttenuationParameters.
        source_parameters (mibitrans.data.read.SourceParameters) : Dataclass object containing source parameters from
            SourceParameters.
        model_parameters (mibitrans.data.read.ModelParameters) : Dataclass object containing model parameters from
            ModelParameters.
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
            biodegradation_capacity (float) : Maximum capacity of contaminant degradation based on electron
            acceptor/donor concentrations, in [g/m^3].
            cxyt_noBC (np.ndarray) : Concentration array in same shape as cxyt, before subtracting
                biodegration_capacity.

        Methods:
            sample : Give concentration at any given position and point in time, closest as discretization allows.

        Raises:
            TypeError : If input is not of the correct Dataclass.

        """
        super().__init__(hydrological_parameters, attenuation_parameters, source_parameters, model_parameters, verbose)
        if auto_run:
            self.cxyt = self._calculate_concentration_for_all_xyt(self.xxx, self.yyy, self.ttt)

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

    def _calculate_biodegradation_capacity(self):
        biodegradation_capacity = 0
        utilization_factor = getattr(self._att_pars, "utilization_factor").dictionary
        for key, item in utilization_factor.items():
            biodegradation_capacity += getattr(self._att_pars, util_to_conc_name[key]) / item

        return biodegradation_capacity

    def _calculate_concentration_for_all_xyt(self, xxx, yyy, ttt):
        cxyt = 0
        with np.errstate(divide="ignore", invalid="ignore"):
            x_term = self._equation_term_x(xxx, ttt)
            additional_x = self._equation_term_additional_x(xxx, ttt)
            z_term = self._equation_term_z(xxx)
            source_decay = self._equation_term_source_decay(xxx, ttt)
            for i in range(len(self.c_source)):
                y_term = self._equation_term_y(i, xxx, yyy)
                cxyt_step = 1 / 8 * self.c_source[i] * source_decay * (x_term + additional_x) * y_term * z_term
                cxyt += cxyt_step

            self.cxyt_noBC = cxyt.copy()
            cxyt -= self.biodegradation_capacity
            cxyt = np.where(cxyt < 0, 0, cxyt)
        self.has_run = True
        return cxyt
