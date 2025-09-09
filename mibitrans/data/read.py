"""Author: Jorrit Bakker.

Module handling data input in the form of a dictionary.
"""
import warnings
import numpy as np
from mibitrans.data.parameter_information import key_dictionary
from dataclasses import dataclass

@dataclass
class HydrologicalParameters:
    """Dataclass handling hydrological parameters."""
    velocity : float = None
    porosity : float = None
    alpha_x : float = None
    alpha_y : float = None
    alpha_z : float = 0
    h_gradient : float = None
    h_conductivity : float = None
    verbose : bool = False

    def __post_init__(self):
        missing_arguments = []
        if self.porosity is None:
            missing_arguments.append("porosity")
        if self.alpha_x is None:
            missing_arguments.append("alpha_x")
        if self.alpha_y is None:
            missing_arguments.append("alpha_y")

        if len(missing_arguments) > 0:
            raise ValueError(f"HydrologicalParameters missing {len(missing_arguments)} arguments: {missing_arguments}.")

        if self.velocity is None and (self.h_gradient is None or self.h_conductivity is None):
            raise ValueError(f"HydrologicalParameters missing required arguments: either velocity or both h_gradient and h_conductivity.")

        # All input parameters should be positive floats, raise appropriate error if not
        for parameter, value in self.__dict__.items():
            # Specific check and error for porosity, which has domain [0,1]
            if parameter == "porosity":
                error = _check_float_fraction(parameter, value)
            elif parameter == "verbose":
                continue
            else:
                error = _check_float_positive(parameter, value)

            if error and (value is not None):
                raise error

        # Velocity is calculated from hydraulic gradient and conductivity when both are given.
        if self.h_gradient and self.h_conductivity:
            # Giving h_gradient and h_conductivity is more specific than giving velocity. So input velocity will be overridden.
            if self.velocity is not None:
                warnings.warn("Both velocity and h_gradient & h_conductivity are defined. Value for velocity will be overridden.", UserWarning)
            self.velocity = self.h_gradient * self.h_conductivity / self.porosity


@dataclass
class AdsorptionDegradationParameters:
    """Dataclass handling adsorption degradation parameters."""
    retardation : float = None
    bulk_density : float = None
    partition_coefficient : float = None
    fraction_organic_carbon : float = None
    decay_rate : float = None
    half_life : float = None
    delta_oxygen : float = None
    delta_nitrate : float = None
    ferrous_iron : float = None
    delta_sulfate : float = None
    methane : float = None
    verbose : bool = False

    def __post_init__(self):
        if self.retardation is None and (self.bulk_density is None
                                         or self.bulk_density is None
                                         or self.partition_coefficient is None
                                         or self.fraction_organic_carbon is None):
            raise ValueError("AdsorptionDegradationParameters missing required arguments: either retardation or (bulk_density, partition_coefficient and fraction_organic_carbon).")

        # All input parameters should be positive floats, raise appropriate error if not
        for parameter, value in self.__dict__.items():
            # Retardation and fraction_organic_carbon have specific domains and are checked separately
            if parameter == "retardation":
                error = _check_float_retardation(parameter, value)
            elif parameter == "fraction_organic_carbon":
                error = _check_float_fraction(parameter, value)
            elif parameter == "verbose":
                continue
            else:
                error = _check_float_positive(parameter, value)

            if error and (value is not None):
                raise error

        # Retardation factor is not calculated from bulk density, partition coefficient and fraction organic carbon
        # in this dataclass, since porosity is required as well, which is defined in the HydrologicalParameters class.

        if self.half_life:
            decay_rate = np.log(2) / self.half_life
            if self.decay_rate and (self.decay_rate != decay_rate):
                warnings.warn("Both contaminant decay rate constant and half life are defined, but are not equal. Only value for decay rate constant will be used in calculations.", UserWarning)
            else:
                self.decay_rate = decay_rate

    def calculate_retardation(self, porosity : float):

        if self.retardation is None:
            self.retardation = 1 + (self.bulk_density / porosity) * self.partition_coefficient * self.fraction_organic_carbon

    def utilization_factor(self):
        """Introduce custom utilization factors for each electron donor/acceptor species."""
        # Come back to later
        return None

@dataclass
class SourceParameters:
    """Dataclass handling source parameters."""
    source_zone_boundary : float | list | np.ndarray = None
    source_zone_concentration : float | list | np.ndarray = None
    depth : float = None
    total_mass : float | str = "infinite"
    verbose : bool = False

    def __post_init__(self):

        # Check if all required arguments are present
        missing_arguments = []
        if self.source_zone_boundary is None:
            missing_arguments.append("source_zone_boundary")
        if self.source_zone_concentration is None:
            missing_arguments.append("source_zone_concentration")
        if self.depth is None:
            missing_arguments.append("depth")

        if len(missing_arguments) > 0:
            raise ValueError(f"SourceParameters missing {len(missing_arguments)} arguments: {missing_arguments}.")

        # Check input value and data type
        for parameter, value in self.__dict__.items():
            if parameter == "verbose":
                continue
            elif parameter == "source_zone_boundary" or parameter == "source_zone_concentration":
                error = _check_array_float_positive(parameter, value)
            # Total source mass can be a positive float or a string denoting that source mass is infinite
            # Therefore, it is checked separately
            elif parameter == "total_mass":
                error = _check_total_mass(parameter, self.total_mass)
            else:
                error = _check_float_positive(parameter, value)

            if error and (value is not None):
                raise error

        # If total mass is any type of string, it is assumed the user's intention is to have it be infinite
        if isinstance(self.total_mass, str):
            self.total_mass = "infinite"

        # Ensure boundary and concentration have same data type
        if isinstance(self.source_zone_boundary, (float, int)):
            self.source_zone_boundary = np.array([self.source_zone_boundary])
        else:
            self.source_zone_boundary = np.array(self.source_zone_boundary)
        if isinstance(self.source_zone_concentration, (float, int)):
            self.source_zone_concentration = np.array([self.source_zone_concentration])
        else:
            self.source_zone_concentration = np.array(self.source_zone_concentration)

        # Each given source zone boundary should have a given concentration, and vice versa
        if self.source_zone_boundary.shape != self.source_zone_concentration.shape:
            raise ValueError(f"Length of source zone boundary ({len(self.source_zone_boundary)}) and source zone concentration ({len(self.source_zone_concentration)}) do not match. Make sure they are of equal length.")

        # Reorder source zone locations if they are not given in order from close to far from source zone center
        if len(self.source_zone_boundary) > 1:
            if not all(self.source_zone_boundary[:-1] <= self.source_zone_boundary[1:]):
                sort_location = np.argsort(self.source_zone_boundary)
                self.source_zone_boundary.sort()
                self.source_zone_concentration = self.source_zone_concentration[sort_location]

                warnings.warn(f"Source zone boundary locations should be ordered by distance from source zone center. Zone boundaries and concentrations have consequently been reordered.")

    def interpolate(self, n_zones, method):
        """Rediscretize source to n zones. Either through linear interpolation or using a normal distribution."""
        # Come back later
        return None


@dataclass
class ModelParameters:
    """Dataclass handling model discretization parameters."""
    model_length : float = None
    model_width : float = None
    model_time : float = None
    dx : float = None
    dy : float = None
    dt : float = None
    verbose : bool = False

    def __post_init__(self):
        # No single argument is required, so no presence check is performed.

        # Check input value and data type
        for parameter, value in self.__dict__.items():
            # Specific check and error for porosity, which has domain [0,1]
            if parameter == "verbose":
                continue
            else:
                error = _check_float_positive(parameter, value)

            if error and (value is not None):
                raise error

        if self.model_length and self.dx:
            if self.model_length / self.dx < 1:
                warnings.warn("Step size is larger than model length, dx will be set to length of model.", UserWarning)
                self.dx = self.model_length
        if self.model_width and self.dy:
            if self.model_width / self.dy < 1:
                warnings.warn("Step size is larger than model width, dy will be set to width of model.", UserWarning)
                self.dy = self.model_width
        if self.model_time and self.dt:
            if self.model_time / self.dt < 1:
                warnings.warn("Step size is larger than model time, dt will be set to time of model.", UserWarning)
                self.dt = self.model_time


def _check_float_positive(parameter : str, value):
    """Check if a variable is a float and if it is positive."""
    if isinstance(value, (float, int)):
        if value >= 0:
            return None
        else:
            return ValueError(f"{parameter} must be >= 0")
    else:
        return TypeError(f"{parameter} must be a float, but is {type(value)} instead.")

def _check_float_fraction(parameter : str, value):
    """Check if a variable is a float and if it is between 0 and 1."""
    if isinstance(value, (float, int)):
        if 0 <= value <= 1:
            return None
        else:
            return ValueError(f"{parameter} must be between 0 and 1")
    else:
        return TypeError(f"{parameter} must be a float, but is {type(value)} instead.")

def _check_float_retardation(parameter : str, value):
    """Check if a variable is a float and if it is 1 or larger."""
    if isinstance(value, (float, int)):
        if value >= 1:
            return None
        else:
            return ValueError(f"{parameter} must be 1 or larger.")
    else:
        return TypeError(f"{parameter} must be a float, but is {type(value)} instead.")

def _check_array_float_positive(parameter : str, value):
    """Check if variable is numpy array, list, or float, if it is positive and if an array is 1-dimensional"""
    if isinstance(value, np.ndarray):
        if len(value.shape) == 1:
            if all(value >= 0):
                return None
            else:
                return ValueError(f"All values in {parameter} should be >= 0.")
        else:
            return ValueError(f"{parameter} should be a float, list or a 1-dimensional array.")

    elif isinstance(value, list):
        if all(isinstance(element, (float, int)) for element in value):
            if all(element >= 0 for element in value):
                return None
            else:
                return ValueError(f"All values in {parameter} should be >= 0.")
        else:
            return TypeError(f"All elements of {parameter} should be a float.")

    elif isinstance(value, (float, int)):
        if value >= 0:
            return None
        else:
            return ValueError(f"{parameter} must be >= 0")

    else:
        return TypeError(f"{parameter} must be a float, list or numpy array, but is {type(value)} instead.")

def _check_total_mass(parameter : str, value):
    """Check variable properties of total source mass specifically."""
    if isinstance(value, (float, int)):
        if value >= 0:
            return None
        else:
            return ValueError(f"{parameter} must be >= 0, or set to 'infinite'.")
    elif isinstance(value, str):
        if value not in ["infinite", "inf", "INF", "Infinite"]:
            warnings.warn(f"{value} is not understood, total source mass is set to infinite.")
        return None
    else:
        return TypeError(f"{parameter} must be a float or 'infinite', but is {type(value)} instead.")

##################
##### Legacy #####
##################

def from_dict(dictionary: dict,
              verbose: bool = True
              ) -> dict:
    """Format and structure input dictionary into a standardized dictionary.

    Args:
        dictionary (dict): Input dictionary.
        verbose (bool, optional): Print verbose output. Defaults to True.

    Returns:
        dict: Dictionary following standardized format.
    """
    params = {}
    unknown_keys = []

    # Convert every input dictionary key by looping over all keys
    for key_input, value_input in dictionary.items():
        key_in_known_keys = False
        for key_params, key_known in key_dictionary.items():
            # Look if input key is listed as a possible name for a parameter
            if key_input in key_known:
                params[key_params] = value_input
                key_in_known_keys = True
                # If input key is recognized, no need to continue this loop and go to next input key.
                break
        if not key_in_known_keys:
            unknown_keys.append(key_input)

    if verbose and len(unknown_keys) > 0:
        print("The following keys were not recognized and not included in output dictionary:", unknown_keys)
    elif verbose:
        print("All keys were recognized")

    return params
