"""Author: Jorrit Bakker.

Module evaluating if a dictionary contains all required (correct) parameters for analysis
"""

import warnings
import numpy as np
import mibitrans.data.parameter_information


def _check_float(parameter: str, value):
    """Check if a variable is a float and if it is positive."""
    if isinstance(value, (float, int, np.floating)):
        return None
    else:
        return TypeError(f"{parameter} must be a float, but is {type(value)} instead.")


def _check_float_positive(parameter: str, value):
    """Check if a variable is a float and if it is positive."""
    is_float = _check_float(parameter, value)
    if is_float is None:
        if value >= 0:
            return None
        else:
            return ValueError(f"{parameter} must be >= 0")
    else:
        return is_float


def _check_float_fraction(parameter: str, value):
    """Check if a variable is a float and if it is between 0 and 1."""
    is_float = _check_float(parameter, value)
    if is_float is None:
        if 0 <= value <= 1:
            return None
        else:
            return ValueError(f"{parameter} must be between 0 and 1")
    else:
        return is_float


def _check_float_retardation(parameter: str, value):
    """Check if a variable is a float and if it is 1 or larger."""
    is_float = _check_float(parameter, value)
    if is_float is None:
        if value >= 1:
            return None
        else:
            return ValueError(f"{parameter} must be 1 or larger.")
    else:
        return TypeError(f"{parameter} must be a float, but is {type(value)} instead.")


def _check_array_float_positive(parameter: str, value):
    """Check if variable is numpy array, list, or float, if it is positive and if an array is 1-dimensional."""
    if isinstance(value, np.ndarray):
        if len(value.shape) == 1:
            if all(value >= 0):
                return None
            else:
                return ValueError(f"All values in {parameter} should be >= 0.")
        else:
            return ValueError(f"{parameter} should be a float, list or a 1-dimensional array.")

    elif isinstance(value, list):
        if all(isinstance(element, (float, int, np.floating)) for element in value):
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


def _check_total_mass(parameter: str, value):
    """Check variable properties of total source mass specifically."""
    if _check_float(parameter, value) is None:
        if value >= 0:
            return None
        else:
            return ValueError(f"{parameter} must be >= 0, or set to 'infinite'.")
    elif isinstance(value, str):
        if value not in ["infinite", "inf", "INF", "Infinite", "Inf"]:
            return ValueError(f"{value} is not understood. For infinite source mass, use 'infinite' or 'inf'.")
        else:
            return None

    else:
        return TypeError(f"{parameter} must be a float or 'infinite', but is {type(value)} instead.")


def _check_dictionary(value):
    """Check if variable is a dictionary, and raise an error if it is not."""
    if not isinstance(value, dict):
        raise TypeError(f"Input must be a dict, but is {type(value)} instead.")


def _check_dataclass(parameter, value, expected_type):
    """Check if variable is of the given type, and raise an error if it is not."""
    if isinstance(value, expected_type):
        return None
    else:
        return TypeError(f"{parameter} must be of type {expected_type}, but is {type(value)} instead.")


def _check_model_type(parameter, allowed_model_types):
    """Check if variable is of the given allowed model types, and raise an error if it is not."""
    if not isinstance(parameter, allowed_model_types):
        raise TypeError(
            f"Input argument model should be in {allowed_model_types.__subclasses__()}, "
            f"but is {type(parameter)} instead."
        )


def _time_check(model, time):
    """Check if time input is valid, and returns the index of nearest time."""
    if time is not None:
        error = _check_float_positive("time", time)
        if error is not None:
            raise error
        elif time > np.max(model.t):
            warnings.warn(
                f"Desired time is larger than maximum time of model ({time} > {np.max(model.t)}). Using maximum time "
                f"of model instead."
            )
            time_pos = -1
        else:
            time_pos = np.argmin(abs(model.t - time))
    else:
        time_pos = -1
    return time_pos


def _y_check(model, y_position):
    """Check if y-position input is valid, and returns the index of nearest y position."""
    error = _check_float("y_position", y_position)
    if error is not None:
        raise error
    if y_position > np.max(model.y):
        warnings.warn(
            f"Desired y position is outside of model domain (abs({y_position}) > {np.max(model.t)}). "
            f"Using closest position inside model domain instead."
        )

    y_pos = np.argmin(abs(model.y - y_position))
    return y_pos


def validate_input_values(parameter, value):
    """Validate if input parameter is of correct type and in correct domain."""
    match parameter:
        # Any input for verbose argument is fine.
        case "verbose":
            error = None
        # Specific check for retardation, which has domain >= 1
        case "retardation":
            error = _check_float_retardation(parameter, value)
        # Specific check for total mass, which can be a positive float, or a specific string
        case "total_mass":
            error = _check_total_mass(parameter, value)
        # Specific check for electron acceptor utilization factor, which should be UtilizationFactor dataclass
        case "utilization_factor":
            error = _check_dataclass(parameter, value, mibitrans.data.parameter_information.UtilizationFactor)
        # Parameters which have domain [0,1]
        case "porosity" | "fraction_organic_carbon":
            error = _check_float_fraction(parameter, value)
        # Parameters which are input as single values, lists or numpy arrays
        case "source_zone_boundary" | "source_zone_concentration":
            error = _check_array_float_positive(parameter, value)
        # All other values are floats on positive domain
        case _:
            error = _check_float_positive(parameter, value)

    if error and (value is not None):
        raise error
