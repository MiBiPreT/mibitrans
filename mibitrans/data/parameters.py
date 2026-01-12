"""Author: Jorrit Bakker.

Module handling data input in the form of a dictionary.
"""

import warnings
from dataclasses import dataclass
import numpy as np
from mibitrans.data.check_input import MissingValueError
from mibitrans.data.check_input import validate_input_values
from mibitrans.data.check_input import validate_source_zones
from mibitrans.visualize.show_conditions import source_zone


@dataclass
class HydrologicalParameters:
    """Dataclass handling input of hydrological parameters.

    Args:
        velocity (float) : Flow velocity in the direction of the groundwater gradient, in [m/d]. Optional if h_gradient
            and h_conductivity are specified.
        h_gradient (float) : Hydraulic gradient of the groundwater, in [m/m]. Optional if velocity is specified.
        h_conductivity (float) : Hydraulic conductivity of the aquifer, in [m/d]. Optional if velocity is specified.
        porosity (float) : Effective soil porosity [-]
        alpha_x (float) : The dispersivity in the x (longitudinal) direction in [m]
        alpha_y (float) : The dispersivity in the y (transverse-horizontal) direction in [m]
        alpha_z (float, optional) : The dispersivity in the z (transverse-vertical) direction in [m]. Defaults to 1e-10
        verbose (bool, optional): Verbose mode. Defaults to False.

    Raises:
        ValueError : If input parameters are incomplete or outside the valid domain.
        TypeError : If input parameters of incorrect datatype.
    """

    velocity: float = None
    h_gradient: float = None
    h_conductivity: float = None
    porosity: float = None
    alpha_x: float = None
    alpha_y: float = None
    alpha_z: float = 1e-10
    verbose: bool = False

    def __setattr__(self, parameter, value):
        """Override parent method to validate input when attribute is set."""
        validate_input_values(parameter, value)
        super().__setattr__(parameter, value)

    def __post_init__(self):
        """Check argument presence, types and domain. Calculate velocity if not given."""
        self._validate_input_presence()

        # Velocity is calculated from hydraulic gradient and conductivity when both are given.
        if self.h_gradient and self.h_conductivity:
            # Giving h_gradient & h_conductivity more specific than giving velocity. Input velocity will be overridden.
            if self.velocity is not None:
                warnings.warn(
                    "Both velocity and h_gradient & h_conductivity are defined. Value for velocity will be overridden.",
                    UserWarning,
                )
            self.velocity = self.h_gradient * self.h_conductivity / self.porosity
            if self.verbose:
                print(f"Groundwater flow velocity has been calculated to be {self.velocity} m/d.")

    def _validate_input_presence(self):
        missing_arguments = []
        if self.porosity is None:
            missing_arguments.append("porosity")
        if self.alpha_x is None:
            missing_arguments.append("alpha_x")
        if self.alpha_y is None:
            missing_arguments.append("alpha_y")

        if len(missing_arguments) > 0:
            raise MissingValueError(
                f"HydrologicalParameters missing {len(missing_arguments)} arguments: {missing_arguments}."
            )

        if self.velocity is None and (self.h_gradient is None or self.h_conductivity is None):
            raise MissingValueError(
                "HydrologicalParameters missing required arguments: either velocity or both h_gradient and"
                "h_conductivity."
            )

        if self.verbose:
            print("All required hydrological input arguments are present.")


@dataclass
class AttenuationParameters:
    """Dataclass handling parameters related to adsorption, diffusion and degradation.

    Args:
        retardation (float) : Retardation factor for transported contaminant [-]. Default is 1.
        decay_rate (float) : First order (linear) decay coefficient in [1/day]. Only required for linear decay models.
            Default is 0. Also sets corresponding half life.
        half_life (float) : Contaminant half life for 1st order (linear) decay, in [days]. Only required for
            linear decay models. Default is 0. Also sets corresponding decay_rate.
        diffusion (float) : Molecular diffusion [m2/day]. Default is 0.
        bulk_density (float) : Soil bulk density, in [g/m^3]. Optional if retardation is specified.
        partition_coefficient (float) : Partition coefficient of the transported contaminant to soil organic matter,
            in [m^3/g]. Optional if retardation is specified.
        fraction_organic_carbon (float) : Fraction of organic material in the soil [-].
            Optional if retardation is specified.
        verbose (bool, optional): Verbose mode. Defaults to False.

    Methods:
        calculate_retardation : Calculate retardation factor from bulk density, partition coefficient and
            fraction organic carbon when given porosity [-]

    Raises:
        ValueError : If input parameters are incomplete or outside the valid domain.
        TypeError : If input parameters of incorrect datatype.
    """

    retardation: float = 1
    decay_rate: float = 0
    half_life: float = 0
    diffusion: float = 0
    bulk_density: float = None
    partition_coefficient: float = None
    fraction_organic_carbon: float = None
    verbose: bool = False

    def __setattr__(self, parameter, value):
        """Override parent method to validate input when attribute is set."""
        validate_input_values(parameter, value)
        # Separate setattr for decay rate and half life because they should always be equivalent
        if parameter == "decay_rate" or parameter == "half_life":
            decay_rate, half_life = self._set_decay(parameter, value)
            super().__setattr__("decay_rate", decay_rate)
            super().__setattr__("half_life", half_life)
        else:
            super().__setattr__(parameter, value)

    def __post_init__(self):
        """Check argument presence, types and domain."""
        self.initialized = True

    def calculate_retardation(self, porosity: float):
        """Calculate retardation factor from soil adsorption parametrers and porosity."""
        self.retardation = (
            1 + (self.bulk_density / porosity) * self.partition_coefficient * self.fraction_organic_carbon
        )
        if self.verbose:
            print(f"Retardation factor has been calculated to be {self.retardation}.")

    def _require_linear_decay(self):
        if self.decay_rate is None and self.half_life is None:
            raise MissingValueError("Linear reaction model requires decay rate or half life.")

    def _set_decay(self, parameter, value):
        if parameter == "decay_rate" and (value != 0 or hasattr(self, "initialized")):
            decay_rate = value
            if value != 0:
                half_life = np.log(2) / value
            else:
                half_life = 0
        elif parameter == "half_life" and (value != 0 or hasattr(self, "initialized")):
            half_life = value
            if value != 0:
                decay_rate = np.log(2) / value
            else:
                decay_rate = 0
        elif parameter == "decay_rate":
            decay_rate = value
            half_life = 0
        elif parameter == "half_life":
            decay_rate = self.decay_rate
            half_life = self.half_life
        else:
            decay_rate = 0
            half_life = 0

        if self.decay_rate != decay_rate and self.decay_rate != 0 and not hasattr(self, "initialized") and value != 0:
            warnings.warn(
                "Both contaminant decay rate and half life were defined, but are not equal. "
                "Value for decay rate will be used.",
                UserWarning,
            )
            half_life = np.log(2) / self.decay_rate
            decay_rate = self.decay_rate

        return decay_rate, half_life


@dataclass
class SourceParameters:
    """Dataclass handling source parameters. Specifying concentrations and extent of source zone.

    Args:
        source_zone_boundary (np.ndarray) : Outer boundary of each source zone, in transverse horizontal direction
            (y-coordiante) [m]. y=0 is at the middle of the contaminant source. Input as numpy array of length equal
            to the amount of source zone. Last value in the array is the limit of the source. For a source with a single
            source zone, only one value is required. Source is symmetrical in the x-axis.
        source_zone_concentration (np.ndarray) : Contaminant concentration in each source zone [g/m^3]. Input as numpy
            array in the same order and of the same length as specified in source_zone_boundary.
        depth (float) : Depth (transverse vertical or z-dimension) of the source zone in [m].
        total_mass (float | str) : Mass of contaminant present in source zone, either expressed in [g],
            or set to 'infinite'. The latter meaning that the source mass and therefore, the source zone concentrations
            do not diminish over time.
        verbose (bool, optional): Verbose mode. Defaults to False.

    Raises:
        ValueError : If input parameters are incomplete or outside the valid domain.
        TypeError : If input parameters of incorrect datatype.
    """

    source_zone_boundary: np.ndarray = None
    source_zone_concentration: np.ndarray = None
    depth: float = None
    total_mass: float | str = "infinite"
    verbose: bool = False

    def __setattr__(self, parameter, value):
        """Override parent method to validate input when attribute is set."""
        validate_input_values(parameter, value)
        if parameter == "total_mass" and (isinstance(value, str) or value == np.inf):
            value = np.inf
        super().__setattr__(parameter, value)
        # When setting source zone boundary or concentration, and both present, check validity in respect to each other.
        if parameter in ["source_zone_boundary", "source_zone_concentration"] and (
            self.source_zone_boundary is not None and self.source_zone_concentration is not None
        ):
            boundary, concentration = validate_source_zones(self.source_zone_boundary, self.source_zone_concentration)
            super().__setattr__("source_zone_boundary", boundary)
            super().__setattr__("source_zone_concentration", concentration)

    def __post_init__(self):
        """Check argument presence, types and domain."""
        self._validate_input_presence()

        # Make sure naming for infinite source mass is consistent from this point onward
        if isinstance(self.total_mass, str):
            self.total_mass = "infinite"

    def interpolate(self, n_zones, method):
        """Rediscretize source to n zones. Either through linear interpolation or using a normal distribution."""
        warnings.warn("This functionality is not implemented yet. Try again later.")
        return None

    def visualize(self):
        """Plot the source zone concentration distribution."""
        source_zone(self)

    def _validate_input_presence(self):
        # Check if all required arguments are present
        missing_arguments = []
        if self.source_zone_boundary is None:
            missing_arguments.append("source_zone_boundary")
        if self.source_zone_concentration is None:
            missing_arguments.append("source_zone_concentration")
        if self.depth is None:
            missing_arguments.append("depth")

        if len(missing_arguments) > 0:
            raise MissingValueError(
                f"SourceParameters missing {len(missing_arguments)} arguments: {missing_arguments}."
            )


@dataclass
class ModelParameters:
    """Dataclass handling model discretization parameters.

    Args:
        model_length (float) : Model extent in the longitudinal (x) direction in [m].
        model_width (float) : Model extent in the transverse horizontal (y) direction in [m].
        model_time (float) : Model duration in [days].
        dx (float, optional) : Model grid discretization step size in the longitudinal (x) direction, in [m]. By
            default, dx = (model_length / 100).
        dy (float, optional) : Model grid discretization step size in the transverse horizontal (y) direction, in [m].
            By default, dy = (model_width / 50).
        dt (float, optional) : Model time discretization step size, in [days]. By default, dt = (model_time / 10).
        verbose (bool, optional): Verbose mode. Defaults to False.

    Raises:
        ValueError : If input parameters are incomplete or outside the valid domain.
        ValueError : If model dimensions are smaller than their given step size.
        TypeError : If input parameters of incorrect datatype.

    """

    model_length: float
    model_width: float
    model_time: float
    dx: float = None
    dy: float = None
    dt: float = None
    verbose: bool = False

    def __setattr__(self, parameter, value):
        """Override parent method to validate input when attribute is set."""
        validate_input_values(parameter, value)
        super().__setattr__(parameter, value)
        self._validate_stepsize(parameter)

    def _validate_stepsize(self, parameter):
        """Validate if model step size is not larger than the corresponding model dimension."""
        match parameter:
            case "dx" | "model_length":
                if self.dx is not None and self.model_length is not None:
                    if self.dx > self.model_length:
                        raise ValueError(
                            f"Model x-direction step size ({self.dx}) "
                            f"is greater than the model length ({self.model_length})."
                        )
            case "dy" | "model_width":
                if self.dy is not None and self.model_width is not None:
                    if self.dy > self.model_width:
                        raise ValueError(
                            f"Model y-direction step size ({self.dy}) "
                            f"is greater than the model width ({self.model_width})."
                        )
            case "dt" | "model_time":
                if self.dt is not None and self.model_time is not None:
                    if self.dt > self.model_time:
                        raise ValueError(
                            f"Model time step size ({self.dt}) "
                            f"is greater than the total model time ({self.model_time})."
                        )

    def __post_init__(self):
        """Set model discretization parameters if not provided."""
        if not self.dx:
            self.dx = self.model_length / 100
        if not self.dy:
            self.dy = self.model_width / 50
        if not self.dt:
            self.dt = self.model_time / 10
