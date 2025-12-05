"""Author: Jorrit Bakker.

Module calculating the mass balance based on base parameters.
"""

import copy
import warnings
import numpy as np
import mibitrans
from mibitrans.analysis.parameter_calculations import calculate_utilization
from mibitrans.data.check_input import check_model_type
from mibitrans.data.check_input import check_time_in_domain


class MassBalance:
    def __init__(self, model, time, method="default", verbose=False):
        self.model = model
        self.verbose = verbose
        self._time_check(time)
        self._method_check(method)

        self._plume_mass_t = None
        self._source_mass_t = None
        self._delta_source_t = None
        self._degraded_mass_t = None
        self._model_without_degradation = None

        # Volume of single cell, as dx * dy * source thickness
        self.cellsize = abs(model.x[0] - model.x[1]) * abs(model.y[0] - model.y[1]) * model.source_parameters.depth

        if self.model.source_parameters.total_mass == "infinite":
            self.source_mass_finite = True
        else:
            self.source_mass_finite = False

        match self.model.mode:
            case "instant_reaction":
                self.model_instant_reaction = True
                self.model_degradation = True
            case "linear" if self.model.attenuation_parameters.decay_rate > 0:
                self.model_degradation = True
                self.model_instant_reaction = False
            case _:
                self.model_degradation = False
                self.model_instant_reaction = False

        self._calculation_routine()

    def __call__(self, time=None, method=None):
        if time:
            self._time_check(time)
        if method:
            self._method_check(method)

        if self.verbose:
            print("Recalculating mass balance...")

        self._calculation_routine()

    @property
    def plume_mass(self):
        return self._plume_mass_t

    @property
    def source_mass(self):
        return self._source_mass_t

    @property
    def delta_source(self):
        return self._delta_source_t

    @property
    def degraded_mass(self):
        return self._degraded_mass_t

    @property
    def model_without_degradation(self):
        return self._model_without_degradation

    def _calculation_routine(self):
        self._plume_mass_t = self.calculate_plume_mass(self.model)
        self._source_mass_t = self.calculate_source_mass()
        self._delta_source_t = self.calculate_delta_source()
        if self.model_degradation:
            self._model_without_degradation = self.calculate_model_without_degradation()
            self._degraded_mass_t = self.calculate_degraded_mass()

    def calculate_plume_mass(self, model):
        # Plume mass of model; concentration is converted to mass by multiplying by cellsize and pore space.
        if isinstance(self.t, np.ndarray):
            _plume_mass_t = np.sum(
                model.cxyt[:, :, 1:] * self.cellsize * self.model.hydrological_parameters.porosity, axis=(1, 2)
            )
        else:
            _plume_mass_t = np.sum(
                model.cxyt[self.t, :, 1:] * self.cellsize * self.model.hydrological_parameters.porosity
            )

        return _plume_mass_t

    def calculate_source_mass(self):
        if isinstance(self.t, np.ndarray):
            local_time = self.t
        else:
            local_time = self.model.t[self.t]

        if self.source_mass_finite:
            _source_mass_t = self.model.source_parameters.source_mass * np.exp(-self.model.k_source * local_time)
        else:
            _source_mass_t = "infinite"

        return _source_mass_t

    def calculate_delta_source(self):
        if isinstance(self.t, np.ndarray):
            local_time = self.t
        else:
            local_time = self.model.t[self.t]

        if self.source_mass_finite:
            _delta_source_t = self.model.source_parameters.source_mass - self._source_mass_t
        else:
            Q, c0_avg = self.model._calculate_discharge_and_average_source_zone_concentration()
            _delta_source_t = Q * c0_avg * local_time
        return _delta_source_t

    def calculate_model_without_degradation(self):
        # Make a no degradation model for comparison, pass the input parameters to a new class instance as kwargs
        _model_without_degradation = self.model.__class__(**self.model.input_parameters)
        _model_without_degradation.attenuation_parameters.decay_rate = 0
        _model_without_degradation.run()
        return _model_without_degradation

    def calculate_degraded_mass(self):
        no_degradation_plume_mass = self.calculate_plume_mass(self._model_without_degradation)
        _degraded_mass_t = no_degradation_plume_mass - self._plume_mass_t
        return _degraded_mass_t

    def _time_check(self, time):
        if time is None or time == "all":
            self.t = self.model.t
        elif isinstance(time, str):
            warnings.warn("String not recognized, defaulting to 'all', for all time points.")
            self.t = self.model.t
        else:
            self.t = check_time_in_domain(self.model, time)

    def _method_check(self, method):
        match method:
            case "default" | 0:
                self.method = "default"
            case "legacy" | 1:
                self.method = "legacy"
            case met if not isinstance(met, (str, int)):
                warnings.warn("Method not recognized, using default.")
                self.method = "default"
            case _:
                raise TypeError(
                    f"Method should a string, either 'default' or 'legacy', but was {type(method)} instead."
                )


def mass_balance(model, time=None) -> dict:
    """Calculate contaminant mass balance across model compartments.

    Args:
        model (mibitrans.transport.model_parent.Transport3D) : Three dimensional transport model object.
        time (float) : Time at which to calculate mass balance. Default is the last time step.

    Returns:
        mass_balance_dict : Dictionary containing the mass balance elements of the given model.
    """
    check_model_type(model, mibitrans.transport.model_parent.Transport3D)
    model = copy.deepcopy(model)
    if not model.has_run:
        model.run()

    time_pos = check_time_in_domain(model, time)

    if isinstance(model.source_parameters.total_mass, str):
        inf_source = True
    else:
        inf_source = False

    mass_balance_dict = {}

    mass_balance_dict["time"] = model.t[time_pos]
    mode = "unknown"
    if hasattr(model, "mode"):
        if model.mode == "instant_reaction":
            mode = model.mode
        elif model.mode == "linear":
            if model._decay_rate == 0:
                mode = "no_decay"
            else:
                mode = "linear_decay"

        if mode in ["instant_reaction", "linear_decay"]:
            no_decay_model = copy.deepcopy(model)
            no_decay_model.mode = "linear"
            no_decay_model.attenuation_parameters.decay_rate = 0
            no_decay_model.run()
        else:
            no_decay_model = model

    # Total source mass at t=0
    M_source_0 = model.source_parameters.total_mass
    mass_balance_dict["source_mass_0"] = M_source_0

    # Total source mass at t=t, for the no decay model
    if inf_source:
        M_source_t = M_source_0
    else:
        M_source_t = M_source_0 * np.exp(-no_decay_model.k_source * model.t[time_pos])
    mass_balance_dict["source_mass_t"] = M_source_t

    # Change in source mass at t=t, due to source decay by transport
    if inf_source:
        Q, c0_avg = no_decay_model._calculate_discharge_and_average_source_zone_concentration()
        M_source_delta = Q * c0_avg * model.t[time_pos]
    else:
        M_source_delta = M_source_0 - M_source_t
    mass_balance_dict["source_mass_change"] = M_source_delta

    # Volume of single cell, as dx * dy * source thickness
    cellsize = abs(model.x[0] - model.x[1]) * abs(model.y[0] - model.y[1]) * model.source_parameters.depth

    # Plume mass of no decay model; concentration is converted to mass by multiplying by cellsize and pore space.
    plume_mass_nodecay = np.sum(
        no_decay_model.cxyt[time_pos, :, 1:] * cellsize * model.hydrological_parameters.porosity
    )
    mass_balance_dict["plume_mass_no_decay"] = plume_mass_nodecay

    # Difference between current plume mass and change in source mass must have been transported outside of model
    # extent for no decay scenarios; preservation of mass.
    if M_source_delta - plume_mass_nodecay < 0:
        transport_outside_extent_nodecay = 0
        mass_balance_dict["transport_outside_extent"] = transport_outside_extent_nodecay
    else:
        transport_outside_extent_nodecay = M_source_delta - plume_mass_nodecay
        mass_balance_dict["transport_outside_extent_nodecay"] = transport_outside_extent_nodecay

    if mode == "linear_decay":
        # Plume mass of linear decay model.
        plume_mass_lindecay = np.sum(model.cxyt[time_pos, :, 1:] * cellsize * model.hydrological_parameters.porosity)
        mass_balance_dict["plume_mass_linear_decay"] = plume_mass_lindecay

        # Calculate transport out of model extent linear decay as fraction of transport out of model for no decay
        # model, scaled by ratio between no decay and linear decay plume mass.
        transport_outside_extent_lindecay = transport_outside_extent_nodecay * plume_mass_lindecay / plume_mass_nodecay
        mass_balance_dict["transport_outside_extent_lineardecay"] = transport_outside_extent_lindecay

        # Contaminant mass degraded by linear decay is difference plume mass no and linear decay plus difference in
        # mass transported outside model extent by no and linear decay.
        degraded_mass = (
            plume_mass_nodecay
            - plume_mass_lindecay
            + transport_outside_extent_nodecay
            - transport_outside_extent_lindecay
        )
        mass_balance_dict["plume_mass_degraded_linear"] = degraded_mass

    elif mode == "instant_reaction":
        # Total source mass at t=t, for the instant reaction model
        if isinstance(model.source_parameters.total_mass, str):
            M_source_t_inst = M_source_0
        else:
            M_source_t_inst = M_source_0 * np.exp(-model.k_source * model.t[time_pos])
        mass_balance_dict["source_mass_instant_t"] = M_source_t_inst

        # Change in source mass at t=t due to source decay by transport and by biodegradation
        if inf_source:
            Q, c0_avg = model._calculate_discharge_and_average_source_zone_concentration()
            M_source_delta = Q * c0_avg * model.t[time_pos]
        else:
            M_source_delta = M_source_0 - M_source_t_inst
        mass_balance_dict["source_mass_instant_change"] = M_source_delta

        # Plume mass without biodegradation according to the instant degradation model
        plume_mass_inst_nodecay = np.sum(
            model.cxyt_noBC[time_pos, :, 1:] * cellsize * model.hydrological_parameters.porosity
        )
        mass_balance_dict["plume_mass_no_decay_instant_reaction"] = plume_mass_inst_nodecay

        # Plume mass with biodegradation according to the instant degradation model
        plume_mass_inst = np.sum(model.cxyt[time_pos, :, 1:] * cellsize * model.hydrological_parameters.porosity)
        mass_balance_dict["plume_mass_instant_reaction"] = plume_mass_inst

        # Assumption: all mass difference between instant degradation model with biodegradation and
        # instant degradation model without biodegradation is caused by degradation.
        degraded_mass = plume_mass_inst_nodecay - plume_mass_inst
        mass_balance_dict["plume_mass_degraded_instant"] = degraded_mass

        # Weight fraction of electron acceptor used for degradation and degraded contaminant
        mass_fraction_electron_acceptor = calculate_utilization(model)

        # Change in total mass of each electron acceptor
        electron_acceptor_mass_change = mass_fraction_electron_acceptor * degraded_mass
        mass_balance_dict["electron_acceptor_mass_change"] = electron_acceptor_mass_change

    return mass_balance_dict
