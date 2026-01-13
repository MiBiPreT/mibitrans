"""Author: Jorrit Bakker.

Module calculating the mass balance based on base parameters.
"""

import warnings
import numpy as np
import mibitrans
from mibitrans.analysis.parameter_calculations import calculate_discharge_and_average_source_zone_concentration
from mibitrans.data.check_input import check_model_type
from mibitrans.data.check_input import check_time_in_domain
from mibitrans.data.check_input import validate_input_values


class MassBalance:
    """Calculate mass balance characteristics of input model."""

    def __init__(self, results, time, verbose=False):
        """Mass balance object with source and plume characteristics at given time(s), of input model.

        Args:
            results: Input model for which mass balance is calculated, should be a child class of Transport3D.
            time (float | str): Time at which to initially calculate the mass balance. Either as a value between 0 and
                model end time. Or as 'all', which will calculate mass balance attributes for each time step as arrays.
            verbose (bool, optional): Verbose mode. Defaults to False.

        Call:
            Calling the MassBalance object will recalculate the mass balance characteristics of input model for given
                input time.

        Properties:
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
        check_model_type(results, mibitrans.transport.model_parent.Results)
        self.results = results
        self.verbose = verbose
        self.t = self._time_check(time)

        self._plume_mass_t = None
        self._source_mass_t = None
        self._delta_source_t = None
        self._degraded_mass_t = None
        self._electron_acceptor_change_t = None
        self._instant_reaction_degraded_mass_t = None
        self._model_without_degradation = None

        # Volume of single cell, as dx * dy * source thickness
        self.cellsize = (
            abs(results.x[0] - results.x[1]) * abs(results.y[0] - results.y[1]) * results.source_parameters.depth
        )

        # Mass balance output differs if the source is represented as an infinite mass.
        if self.results.source_parameters.total_mass == np.inf:
            self.source_mass_finite = False
        else:
            self.source_mass_finite = True

        match self.results.mode:
            # Instant reaction model
            case "instant_reaction":
                self.model_instant_reaction = True
                self.model_degradation = True
            # Linear decay model
            case "linear" if self.results.attenuation_parameters.decay_rate > 0:
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
            self.t = self._time_check(time)

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
            time_to_threshold = (
                -1 / self.results.k_source * np.log(threshold / self.results.source_parameters.total_mass)
            )
        return time_to_threshold

    def _calculation_routine(self):
        """Perform mass_balance calculations."""
        self._check_model_extent()
        self._plume_mass_t = self._calculate_plume_mass(self.results)
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
        # Relative concentration considered to be boundary of the plume extent.
        extent_threshold_value = 0.01
        if isinstance(self.t, np.ndarray):
            cxyt_y_boundary = self.results.relative_cxyt[:, [0, -1], :]
            cxyt_x_boundary = self.results.relative_cxyt[:, :, -1]
        else:
            cxyt_y_boundary = self.results.relative_cxyt[self._t_index, [0, -1], :]
            cxyt_x_boundary = self.results.relative_cxyt[self._t_index, :, -1]

        y_boundary_above_threshold = np.where(cxyt_y_boundary > extent_threshold_value, cxyt_y_boundary, 0.0)
        x_boundary_above_threshold = np.where(cxyt_x_boundary > extent_threshold_value, cxyt_x_boundary, 0.0)
        if np.sum(y_boundary_above_threshold) > 0:
            y_max = np.round(
                np.max(y_boundary_above_threshold) * np.max(self.results.source_parameters.source_zone_concentration), 2
            )
            warnings.warn(
                "Contaminant plume extents beyond the model width, with a maximum concentration at the "
                f"boundary of {y_max}g/m3. To ensure reliable mass balance, re-run the model with increased dimensions "
                "to include the entire plume width in the model extent."
            )
        if np.sum(x_boundary_above_threshold) > 0:
            x_max = np.round(
                np.max(x_boundary_above_threshold) * np.max(self.results.source_parameters.source_zone_concentration), 2
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
            plume_mass_t = np.sum(
                model.cxyt[:, :, 1:] * self.cellsize * self.results.hydrological_parameters.porosity, axis=(1, 2)
            )
        else:
            plume_mass_t = np.sum(
                model.cxyt[self._t_index, :, 1:] * self.cellsize * self.results.hydrological_parameters.porosity
            )

        return plume_mass_t

    def _calculate_source_mass(self):
        """Calculate source mass of input model, for the given time(s)."""
        if self.source_mass_finite:
            source_mass_t = self.results.source_parameters.total_mass * np.exp(-self.results.k_source * self.t)
        else:
            source_mass_t = np.inf

        return source_mass_t

    def _calculate_delta_source(self):
        """Calculate difference in source mass between t=0 and given time(s)."""
        if self.source_mass_finite:
            delta_source_t = self.results.source_parameters.total_mass - self._source_mass_t
        else:
            Q, c0_avg = calculate_discharge_and_average_source_zone_concentration(self.results)
            delta_source_t = Q * c0_avg * self.t
        return delta_source_t

    def _calculate_model_without_degradation(self):
        """Make a no degradation model for comparison, pass the input parameters to a new class instance as kwargs."""
        model_without_degradation = self.results.model_type(**self.results.input_parameters)
        model_without_degradation.attenuation_parameters.decay_rate = 0
        model_without_degradation.run()
        return model_without_degradation

    def _calculate_degraded_mass(self):
        """Calculate difference between input model plume mass and no degradation model, for a given time(s)."""
        no_degradation_plume_mass = self._calculate_plume_mass(self._model_without_degradation)
        degraded_mass_t = no_degradation_plume_mass - self._plume_mass_t
        return degraded_mass_t

    def _calculate_instant_reaction_degraded_mass(self):
        """Calculate difference between input model plume mass and no degradation model for instant reaction model."""
        if isinstance(self.t, np.ndarray):
            plume_mass_t_noBC = np.sum(
                self.results.cxyt_noBC[:, :, 1:] * self.cellsize * self.results.hydrological_parameters.porosity,
                axis=(1, 2),
            )
        else:
            plume_mass_t_noBC = np.sum(
                self.results.cxyt_noBC[self._t_index, :, 1:]
                * self.cellsize
                * self.results.hydrological_parameters.porosity
            )

        degraded_mass_instant_t = plume_mass_t_noBC - self._plume_mass_t
        return degraded_mass_instant_t

    def _calculate_electron_acceptor_change(self):
        """Calculate the change in electron acceptor mass for given time(s) for instant reaction model."""
        mass_fraction_degraded_acceptor = self.results.electron_acceptors.array / self.results.biodegradation_capacity
        electron_acceptor_change = {}
        electron_acceptors = ["oxygen", "nitrate", "ferrous_iron", "sulfate", "methane"]
        for i, ea in enumerate(electron_acceptors):
            electron_acceptor_change[ea] = self._instant_reaction_degraded_mass_t * mass_fraction_degraded_acceptor[i]

        return electron_acceptor_change

    def _time_check(self, time):
        """Check if time input is valid."""
        if time is None or time == "all":
            t = self.results.t
        elif isinstance(time, str):
            warnings.warn("String not recognized, defaulting to 'all', for all time points.")
            t = self.results.t
        else:
            self._t_index = check_time_in_domain(self.results, time)
            t = float(self.results.t[self._t_index])
        return t
