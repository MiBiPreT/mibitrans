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
        print(self.__dict__)

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

    def utilization_factor(self):
        """Introduce custom utilization factors for each electron donor/acceptor species."""
        # Come back to later
        return None

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

        for parameter, value in self.__dict__.items():
            # Specific check and error for porosity, which has domain [0,1]
            if parameter == "verbose":
                continue
            else:
                error = _check_float_positive(parameter, value)

            if error and (value is not None):
                raise error




class SourceParameters:
    """Dataclass handling source parameters."""
    source_zone_boundary : float | np.ndarray = None
    source_zone_concentration : float | np.ndarray  = None
    depth : float = None
    total_mass : float | str = None

    def __post_init__(self):
        return None

    def interpolate(self, n_zones, method):
        """Rediscretize source to n zones. Either through linear interpolation or using a normal distribution."""
        return None

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


