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
from mibitrans.data.check_input import validate_input_values


class MassBalance:
    """Calculate mass balance characteristics of input model."""

    def __init__(self, model, time, method, verbose=False):
        """Mass balance object with source and plume characteristics at given time(s), of input model.

        Args:
            model: Input model for which mass balance is calculated, should be a child class of Transport3D.
            time (float | str): Time at which to initially calculate the mass balance. Either as a value between 0 and
                model end time. Or as 'all', which will calculate mass balance attributes for each time step as arrays.
            method (str): Which method to use for mass balance. 'default' will generate a mass balance based on model
                in and output. 'legacy' will generate a mass balance as done in older mibitrans versions, which is based
                on BIOSCREEN mass balance. Note that the latter is not conservative in inferences made on data, and
                should be used with discretion.
            verbose (bool, optional): Verbose mode. Defaults to False.

        Call:
            Calling the MassBalance object will recalculate the mass balance characteristics of input model for given
                input time and method.

        Properties:
            plume_mass: Mass of the contaminant plume inside the model extent, at the given time(s), in [g].
            source_mass: Mass of the contaminant source at the given time(s), in [g]. No values are given for models
                with infinite source mass.
            delta_source: Difference in mass between contaminant source at given time and source at t = 0, in [g].
            degraded_mass: Mass of plume contaminant degradation at the given time(s), compared to a model without
                degradation, in [g]. Has no value if model does not consider degradation.
            model_without_degradation: Object of model without degradation. Has no value if model does not consider
                degradation.
        """
        check_model_type(model, mibitrans.transport.model_parent.Transport3D)
        self.model = model
        self.verbose = verbose
        self._time_check(time)
        self._method_check(method)

        self._plume_mass_t = None
        self._source_mass_t = None
        self._delta_source_t = None
        self._degraded_mass_t = None
        self._electron_acceptor_change_t = None
        self._instant_reaction_degraded_mass_t = None
        self._model_without_degradation = None

        # Relative concentration considered to be boundary of the plume extent.
        self.extent_threshold_value = 0.01
        # Volume of single cell, as dx * dy * source thickness
        self.cellsize = abs(model.x[0] - model.x[1]) * abs(model.y[0] - model.y[1]) * model.source_parameters.depth

        # Mass balance output differs if the source as infinite mass
        if self.model.source_parameters.total_mass == "infinite":
            self.source_mass_finite = False
        else:
            self.source_mass_finite = True

        match self.model.mode:
            # Instant reaction model
            case "instant_reaction":
                self.model_instant_reaction = True
                self.model_degradation = True
            # Linear decay model
            case "linear" if self.model.attenuation_parameters.decay_rate > 0:
                self.model_degradation = True
                self.model_instant_reaction = False
            # No decay model
            case _:
                self.model_degradation = False
                self.model_instant_reaction = False

        if self.verbose:
            print("Calculating mass balance...")

        self._calculation_routine()

    def __call__(self, time=None, method=None):
        """Recalculate the mass balance characteristics of input model for given time and method."""
        if time:
            self._time_check(time)
        if method:
            self._method_check(method)

        if self.verbose:
            print("Recalculating mass balance...")

        self._calculation_routine()

    @property
    def plume_mass(self):
        """Mass of the contaminant plume in the model extent, at the given time(s), in [g]."""
        return self._plume_mass_t

    @property
    def source_mass(self):
        """Mass of the contaminant source at the given time(s), in [g]. No values are given for infinite source mass."""
        return self._source_mass_t

    @property
    def delta_source(self):
        """Difference in mass between contaminant source at given time and source at t = 0, in [g]."""
        return self._delta_source_t

    @property
    def degraded_mass(self):
        """Mass of plume contaminant degradation at the given time(s), compared to a no degradation model, in [g]."""
        return self._degraded_mass_t

    @property
    def model_without_degradation(self):
        """Model with no degradation used to compare with given model."""
        return self._model_without_degradation

    @property
    def instant_reaction_degraded_mass(self):
        """Difference in plume mass instant reaction with and without biodegradation capacity subtracted, in [g].

        For the instant reaction model, the underlying assumption reads that observed concentrations in the source zone
        are post-degradation. Therefore, the source concentrations without any biodegradation would be higher, the
        amount which is determined by the biodegradation capacity. Then, according to this method, the degraded mass
        is the difference between plume mass before and after subtracting the biodegradation capacity.
        """
        return self._instant_reaction_degraded_mass_t

    @property
    def electron_acceptor_change(self):
        """Change in electron acceptor/byproduct masses at the given time(s), in [g]. Only for instant reaction.

        Electron acceptor/byproduct consumption or generation is based on the degraded plume mass (specifically
        'instant_reaction_degraded_mass'), the utilization factor and relative abundance of the acceptors/byproducts.
        Under the governing assumptions of the instant reaction model, a crude estimate of the total consumption of
        electron acceptors and the generation of byproduct is calculated.
        """
        return self._electron_acceptor_change_t

    def source_threshold(self, threshold):
        """Calculate when source mass is below given threshold. No values are given for infinite source mass."""
        validate_input_values("threshold", threshold)
        if not self.source_mass_finite:
            raise ValueError("Source mass is infinite and therefore cannot go below given threshold.")
        else:
            time_to_threshold = -1 / self.model.k_source * np.log(threshold / self.model.source_parameters.total_mass)
        return time_to_threshold

    def _calculation_routine(self):
        """Performs mass_balance calculations"""
        self._check_model_extent()
        self._plume_mass_t = self._calculate_plume_mass(self.model)
        self._source_mass_t = self._calculate_source_mass()
        self._delta_source_t = self._calculate_delta_source()
        if self.model_degradation:
            self._model_without_degradation = self._calculate_model_without_degradation()
            self._degraded_mass_t = self._calculate_degraded_mass()
        if self.model_instant_reaction:
            self._instant_reaction_degraded_mass_t = self._calculate_instant_reaction_degraded_mass()
            self._electron_acceptor_change_t = self._calculate_electron_acceptor_change()

    def _check_model_extent(self):
        """Check if contaminant plume at given time is reasonably situated within the model extent."""
        if isinstance(self.t, np.ndarray):
            cxyt_y_boundary = self.model.relative_cxyt[:, [0, -1], :]
            cxyt_x_boundary = self.model.relative_cxyt[:, :, -1]
        else:
            cxyt_y_boundary = self.model.relative_cxyt[self.t, [0, -1], :]
            cxyt_x_boundary = self.model.relative_cxyt[self.t, :, -1]

        y_boundary_above_threshold = np.where(cxyt_y_boundary > self.extent_threshold_value, cxyt_y_boundary, 0.0)
        x_boundary_above_threshold = np.where(cxyt_x_boundary > self.extent_threshold_value, cxyt_x_boundary, 0.0)
        if np.sum(y_boundary_above_threshold) > 0:
            y_max = np.round(
                np.max(y_boundary_above_threshold) * np.max(self.model.source_parameters.source_zone_concentration), 2
            )
            warnings.warn(
                "Contaminant plume extents beyond the model width, with a maximum concentration at the "
                f"boundary of {y_max}g/m3. To ensure reliable mass balance, re-run the model with increased dimensions "
                "to include the entire plume width in the model extent."
            )
        if np.sum(x_boundary_above_threshold) > 0:
            x_max = np.round(
                np.max(x_boundary_above_threshold) * np.max(self.model.source_parameters.source_zone_concentration), 2
            )
            warnings.warn(
                "Contaminant plume extents beyond the model length, with a maximum concentration at the "
                f"boundary of {x_max}g/m3. To ensure reliable mass balance, re-run the model with increased dimensions "
                "to include the entire plume length in the model extent."
            )

    def _calculate_plume_mass(self, model):
        """Calculate plume mass of input model, for the given time(s)."""
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

    def _calculate_source_mass(self):
        """Calculate source mass of input model, for the given time(s)."""
        if isinstance(self.t, np.ndarray):
            local_time = self.t
        else:
            local_time = self.model.t[self.t]

        if self.source_mass_finite:
            _source_mass_t = self.model.source_parameters.total_mass * np.exp(-self.model.k_source * local_time)
        else:
            _source_mass_t = "infinite"

        return _source_mass_t

    def _calculate_delta_source(self):
        """Calculate difference in source mass between t=0 and given time(s)."""
        if isinstance(self.t, np.ndarray):
            local_time = self.t
        else:
            local_time = self.model.t[self.t]

        if self.source_mass_finite:
            _delta_source_t = self.model.source_parameters.total_mass - self._source_mass_t
        else:
            Q, c0_avg = self.model._calculate_discharge_and_average_source_zone_concentration()
            _delta_source_t = Q * c0_avg * local_time
        return _delta_source_t

    def _calculate_model_without_degradation(self):
        """Make a no degradation model for comparison, pass the input parameters to a new class instance as kwargs."""
        _model_without_degradation = self.model.__class__(**self.model.input_parameters)
        _model_without_degradation.attenuation_parameters.decay_rate = 0
        _model_without_degradation.run()
        return _model_without_degradation

    def _calculate_degraded_mass(self):
        """Calculate difference between input model plume mass and no degradation model, for a given time(s)."""
        no_degradation_plume_mass = self._calculate_plume_mass(self._model_without_degradation)
        _degraded_mass_t = no_degradation_plume_mass - self._plume_mass_t
        return _degraded_mass_t

    def _calculate_instant_reaction_degraded_mass(self):
        if isinstance(self.t, np.ndarray):
            _plume_mass_t_noBC = np.sum(
                self.model.cxyt_noBC[:, :, 1:] * self.cellsize * self.model.hydrological_parameters.porosity,
                axis=(1, 2),
            )
        else:
            _plume_mass_t_noBC = np.sum(
                self.model.cxyt_noBC[self.t, :, 1:] * self.cellsize * self.model.hydrological_parameters.porosity
            )

        degraded_mass_instant_t = _plume_mass_t_noBC - self._plume_mass_t
        return degraded_mass_instant_t

    def _calculate_electron_acceptor_change(self):
        mass_fraction_degraded_acceptor = self.model._electron_acceptors.array / self.model.biodegradation_capacity
        electron_acceptor_change = {}
        electron_acceptors = ["oxygen", "nitrate", "ferrous_iron", "sulfate", "methane"]
        for i, ea in enumerate(electron_acceptors):
            electron_acceptor_change[ea] = self._instant_reaction_degraded_mass_t * mass_fraction_degraded_acceptor[i]

        return electron_acceptor_change

    def _time_check(self, time):
        """Check if time input is valid."""
        if time is None or time == "all":
            self.t = self.model.t
        elif isinstance(time, str):
            warnings.warn("String not recognized, defaulting to 'all', for all time points.")
            self.t = self.model.t
        else:
            self.t = check_time_in_domain(self.model, time)

    def _method_check(self, method):
        """Check if method input is valid."""
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
