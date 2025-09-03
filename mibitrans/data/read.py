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

        for parameter, value in self.__dict__.items():
            error = _check_positive_float(parameter, value)
            if error and (value is not None):
                raise error
        if self.porosity > 1:
            raise ValueError(f"HydrologicalParameters porosity must be between 0 and 1.")

        if self.h_gradient and self.h_conductivity:
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

class ModelParameters:
    """Dataclass handling model discretization parameters."""
    model_length : float = None
    model_width : float = None
    model_time : float = None
    dx : float = None
    dy : float = None
    dt : float = None

class SourceParameters:
    """Dataclass handling source parameters."""
    source_zone_boundary : float | np.ndarray = None
    source_zone_concentration : float | np.ndarray  = None
    depth : float = None
    total_mass : float = None

    def interpolate(self, n_zones, method):
        """Rediscretize source to n zones."""
        return None

def _check_positive_float(parameter : str, value):
    """Check if a variable is a float and if it is positive."""
    if isinstance(value, (float, int)):
        if value >= 0:
            return None
        else:
            return ValueError(f"{parameter} must be >= 0")
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


