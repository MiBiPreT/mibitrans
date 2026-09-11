"""Author: Jorrit Bakker.

Module evaluating if a dictionary contains all required (correct) parameters for analysis
"""

import warnings
import numpy as np
import mibitrans


# General input checking function
def validate_input_values(parameter, value, expectation=None):
    """Validate if input parameter is of correct type and in correct domain."""
    match parameter:
        # Any input for verbose argument is fine; if set to anything other than False, 0 or None, verbose is on
        case "verbose":
            error = None
        case "_on_change":
            error = None
        # Specific check for retardation, which has domain >= 1
        case "retardation":
            error = _check_numeric_retardation(parameter, value)
        # Specific check for total mass, which can be a positive float, or a specific string
        case "total_mass":
            error = _check_total_mass(parameter, value)
        # Specific check for electron acceptor utilization factor, which should be UtilizationFactor dataclass
        case "utilization_factor":
            error = _check_dataclass(parameter, value, mibitrans.data.parameter_information.UtilizationFactor)
        case "hydrological_parameters" | "attenuation_parameters" | "source_parameters" | "model_parameters":
            error = _check_dataclass(parameter, value, expectation)
        # Parameters which can be any float value
        case "y_position":
            error = _check_numeric(parameter, value)
        # Parameters which have domain [0,1]
        case "porosity" | "fraction_organic_carbon":
            error = _check_numeric_fraction(parameter, value)
        # Parameters which are input as single values, lists or numpy arrays
        case "source_zone_boundary" | "decay_rate" | "half_life" | "mass_ratios":
            error = _check_array_list_numeric_positive(parameter, value, sublist_allowed=False)
        # Parameters which are input as single values, lists or numpy arrays, and may contain nested lists/arrays
        case "source_zone_concentration":
            error = _check_array_list_numeric_positive(parameter, value, sublist_allowed=True)
        case "electron_acceptors":
            error = _check_electron_acceptor(value)
        # All other parameters are checked as floats on positive domain
        case _:
            error = _check_numeric_positive(parameter, value)

    if error and (value is not None):
        raise error


# Protected checking functions used by `validate_input_values`
def _check_numeric(parameter: str, value):
    """Check if a variable is numerical and if it is positive."""
    if isinstance(value, (float, int, np.floating, np.integer)):
        return None
    else:
        return TypeError(f"{parameter} must be a float, but is {type(value)} instead.")


def _check_numeric_positive(parameter: str, value):
    """Check if a variable is numerical and if it is positive."""
    is_float = _check_numeric(parameter, value)
    if is_float is None:
        if value >= 0:
            return None
        else:
            return DomainValueError(parameter, ">=0", value)
    else:
        return is_float


def _check_numeric_fraction(parameter: str, value):
    """Check if a variable is numerical and if it is between 0 and 1."""
    is_float = _check_numeric(parameter, value)
    if is_float is None:
        if 0 <= value <= 1:
            return None
        else:
            return DomainValueError(parameter, "0<n<1", value)
    else:
        return is_float


def _check_numeric_retardation(parameter: str, value):
    """Check if a variable is numerical and if it is 1 or larger."""
    is_float = _check_numeric(parameter, value)
    if is_float is None:
        if value >= 1:
            return None
        else:
            return DomainValueError(parameter, ">=1", value)
    else:
        return TypeError(f"{parameter} must be a float, but is {type(value)} instead.")


def _check_array_list_numeric_positive(parameter: str, value, sublist_allowed):
    """Check if variable is numpy array, list, or numerical, if it is positive and if an array is 1-dimensional."""
    if isinstance(value, np.ndarray) and len(value.shape) != 1:
        return ValueError(
            f"{parameter} must be a 1D array/list of floats or list of 1D-arrays/floats, not a multi-dimensional array."
        )
    if isinstance(value, (np.ndarray, list)) and sublist_allowed:
        return _check_nested_list_array_positive(parameter, value)
    elif isinstance(value, (np.ndarray, list)):
        return _check_list_array_positive(parameter, value, sublist_allowed)
    elif isinstance(value, (float, int, np.floating, np.integer)):
        return _check_numeric_positive(parameter, value)
    else:
        return TypeError(f"{parameter} must be a float, list or numpy array, but is {type(value)} instead.")


def _check_nested_list_array_positive(parameter: str, value):
    """Check if nested arrays/lists only contain positive numeric values."""
    if not all(isinstance(item, (list, np.ndarray)) for item in value):
        return _check_list_array_positive(parameter, value, True)
    for item in value:
        return_value = _check_list_array_positive(parameter, item)
        if return_value:
            return return_value
    return None


def _check_list_array_positive(parameter: str, value, sublist_allowed=False):
    """Check if arrays/lists only contain positive numeric values."""
    if not all(isinstance(item, (int, float, np.floating, np.integer)) for item in value):
        error_message = f"All (sub-)elements of {parameter} should be a float."
        if sublist_allowed:
            error_message += " Or a list/array."
        return TypeError(error_message)
    if not all(item >= 0 for item in value):
        return DomainValueError(parameter, ">=0", value, prefix="All (sub-)elements in")
    else:
        return None


def validate_source_zones(boundary, concentration):
    """Validate and adapt input of source_zone_boundary and source_zone_concentration arrays."""
    # To make sure that broadcasting and indexing works properly in the models, transform everything into arrays.
    # Ensure source_zone_boundary is a numpy array
    boundary = _check_source_boundary_as_array(boundary)
    # Ensure that source_zone_concentration is a numpy array, or a list of numpy arrays for chain decay.
    if isinstance(concentration, list):
        concentration = _check_source_concentrations_as_arrays(concentration)
    elif isinstance(concentration, (int, float, np.floating, np.integer)):
        concentration = np.array([concentration])

    # When only a single source zone boundary is given, but multiple source zone concentrations, it could be
    # interpreted as invalid input for a source zone of a single contaminant. However, it could also be considered as
    # multiple single source concentrations for multiple contaminants in chain decay. Therefore, no
    # error will be raised if length of boundary array != length concentration array.
    if len(boundary) != len(concentration) and len(boundary) == 1:
        concentration = [np.array([conc]) for conc in concentration]

    # Source zone boundary should be ordered by distance from source center to fringes, to make source zone input
    # less ambiguous
    _check_source_boundary_order(boundary)

    check_conc = [concentration] if not isinstance(concentration, list) else concentration
    for conc in check_conc:
        # Each given source zone boundary should have a corresponding concentration, and vice versa, of equal length
        _check_source_boundary_concentration_length(boundary, conc)
        # Superposition method as implemented only works if a zone closer to the center has higher concentration than
        # outer zones.
        _check_source_concentration_order(conc)
    return boundary, concentration


def _check_source_boundary_as_array(boundary):
    """Ensure that source zone boundary is of the type np.ndarray."""
    if isinstance(boundary, (float, int, np.floating, np.integer)):
        return np.array([boundary])
    else:
        return np.array(boundary)


def _check_source_concentrations_as_arrays(concentration):
    """Ensure that source zone concentration is of the type np.ndarray or list(np.ndarray)."""
    if isinstance(concentration[0], (list, np.ndarray)):
        if len(concentration) == 1:
            concentration = np.array(concentration[0])
        else:
            concentration = [np.array(conc) for conc in concentration]
    else:
        concentration = np.array(concentration)
    return concentration


def _check_source_boundary_order(boundary: np.ndarray):
    """Check if source zone boundary is ordered from low to high values."""
    if len(boundary) > 1 and (not all(boundary[:-1] <= boundary[1:])):
        boundary.sort()
        raise ValueError(
            "source_zone_boundary locations should be ordered by distance from source zone center. Thus, current "
            f"input for source_zone_boundary is supposed to be {boundary}. source_zone_concentration should be re-"
            f"ordered accordingly as well, with highest concentrations at the innermost source zone."
        )


def _check_source_boundary_concentration_length(boundary: np.ndarray, conc: np.ndarray):
    """Check if amount of source boundaries corresponds with amount of source zone concentrations."""
    if boundary.shape != conc.shape:
        try:
            len_conc = len(conc)
        except TypeError:
            conc = np.array([conc])
            len_conc = len(conc)

        raise ValueError(
            f"Length of source zone boundary (len={len(boundary)}, for {boundary}) and source zone concentration "
            f"(len={len_conc}, for {conc}) do not match. Make sure they are of equal length."
        )


def _check_source_concentration_order(conc):
    """Check if source zone concentration is ordered from high to low."""
    if not all(conc[:-1] >= conc[1:]) and not isinstance(conc, (float, int, np.floating, np.integer)):
        raise ValueError(
            "Source zone concentrations should be in descending order; no source zone can have a concentration "
            "higher than the concentration of a zone closer to source center, due to the superposition method."
        )


def _check_total_mass(parameter: str, value):
    """Check variable properties of total source mass specifically."""
    if isinstance(value, str):
        if "inf" not in value:
            return ValueError(f"{value} is not understood. For infinite source mass, use 'infinite' or 'np.inf'.")
        else:
            return None
    elif _check_numeric(parameter, value) is None:
        if value >= 0:
            return None
        else:
            return DomainValueError(parameter, ">=0, 'infinite', or np.inf", value)
    else:
        return TypeError(f"{parameter} must be a float or 'infinite', but is {type(value)} instead.")


def _check_dataclass(parameter, value, expected_type):
    """Check if variable is of the given type, and raise an error if it is not."""
    if isinstance(value, expected_type):
        return None
    else:
        return TypeError(f"{parameter} must be of type {expected_type}, but is {type(value)} instead.")


def _check_electron_acceptor(value):
    """Check if variable is an ElectronAcceptors dataclass, list, array or dictionary, raise an error if it is not."""
    if isinstance(value, mibitrans.data.parameter_information.ElectronAcceptors):
        return None
    elif isinstance(value, (list, np.ndarray, dict)):
        if len(value) != 5:
            return ValueError(
                f"Input for electron_acceptors as list, array or dictionary must have an entry for each electron "
                f"acceptor, of which there are five utilized by this model. The current input has {len(value)} "
                f"entries instead."
            )
        else:
            return None
    else:
        return TypeError(
            f"electron_acceptors must be of type {mibitrans.data.parameter_information.ElectronAcceptors},"
            f" or alternatively as list, numpy array or dictionary containing electron acceptor "
            f"concentrations. But is {type(value)} instead."
        )


# Unprotected checking functions


def check_instant_reaction_acceptor_input(electron_acceptors, utilization_factor):
    """Check if electron acceptor and utilization factor are of correct datatype. Then pass them to dataclasses."""
    if isinstance(electron_acceptors, (list, np.ndarray)):
        electron_acceptors_out = mibitrans.data.parameter_information.ElectronAcceptors(*electron_acceptors)
    elif isinstance(electron_acceptors, dict):
        electron_acceptors_out = mibitrans.data.parameter_information.ElectronAcceptors(**electron_acceptors)
    elif isinstance(electron_acceptors, mibitrans.data.parameter_information.ElectronAcceptors):
        electron_acceptors_out = electron_acceptors
    else:
        raise TypeError(
            f"electron_acceptors must be a list, dictionary or ElectronAcceptors dataclass, but is "
            f"{type(electron_acceptors)} instead."
        )

    if isinstance(utilization_factor, (list, np.ndarray)):
        utilization_factor_out = mibitrans.data.parameter_information.UtilizationFactor(*utilization_factor)
    elif isinstance(utilization_factor, dict):
        utilization_factor_out = mibitrans.data.parameter_information.UtilizationFactor(**utilization_factor)
    elif isinstance(utilization_factor, mibitrans.data.parameter_information.UtilizationFactor):
        utilization_factor_out = utilization_factor
    else:
        raise TypeError(
            f"utilization_factor must be a list, dictionary or UtilizationFactor dataclass, but is "
            f"{type(utilization_factor)} instead."
        )

    return electron_acceptors_out, utilization_factor_out


def check_chain_decay_validity(attenuation_parameters, source_parameters, mass_ratios):
    """Check if input for chain decay is valid."""
    if not attenuation_parameters.chain_decay:
        raise ValueError(
            "Attenuation parameters does not contain information for chain decay. Decay rate should be "
            "provided as list or array of degradation rates."
        )
    if not source_parameters.chain_decay_source:
        raise ValueError(
            "Source parameters does not contain information for chain decay. Separate source zone "
            "concentrations should be given for each compound in the chain."
        )
    if len(mass_ratios) != len(attenuation_parameters.decay_rate) - 1:
        raise ValueError(
            "Length of mass_ratios array should be one less than length of decay_rate. As for the "
            "degradation of the final compound in the chain decay, mass ratio is irrelevant."
        )
    if len(source_parameters.source_zone_concentration) != len(attenuation_parameters.decay_rate):
        raise ValueError(
            "Amount of provided source zone concentrations should be equal to the amount of provided decay rates."
            "One for each compound in the chain decay."
        )


def check_model_type(parameter, allowed_model_types):
    """Check if variable is of the given allowed model types, and raise an error if it is not."""
    if not isinstance(parameter, allowed_model_types):
        if isinstance(allowed_model_types, tuple):
            raise TypeError(
                f"Input argument model should be subclass of {allowed_model_types}, but is {type(parameter)} instead."
            )
        else:
            raise TypeError(
                f"Input argument model should be in {allowed_model_types.__subclasses__()}, "
                f"but is {type(parameter)} instead."
            )


def check_x_in_domain(model, x_position):
    """Check if x-position input is valid, and returns the index of nearest x position."""
    error = _check_numeric_positive("x_position", x_position)
    if error is not None:
        raise error
    if x_position > np.max(model.x):
        warnings.warn(
            f"Desired x position is outside of model domain ({x_position} > {np.max(model.x)}). "
            f"Using closest position inside model domain instead."
        )

    x_pos = np.argmin(abs(model.x - x_position))
    return x_pos


def check_y_in_domain(model, y_position):
    """Check if y-position input is valid, and returns the index of nearest y position."""
    error = _check_numeric("y_position", y_position)
    if error is not None:
        raise error
    if y_position > np.max(model.y):
        warnings.warn(
            f"Desired y position is outside of model domain (abs({y_position}) > {np.max(model.y)}). "
            f"Using closest position inside model domain instead."
        )

    y_pos = np.argmin(abs(model.y - y_position))
    return y_pos


def check_time_in_domain(model, time):
    """Check if time input is valid, and returns the index of nearest time."""
    if time is not None:
        error = _check_numeric_positive("time", time)
        if error is not None:
            raise error
        elif time > np.max(model.t):
            warnings.warn(
                f"Desired time is larger than maximum time of model ({time} > {np.max(model.t)}). Using maximum time "
                f"of model instead."
            )
            time_pos = len(model.t) - 1
        else:
            time_pos = np.argmin(abs(model.t - time))
    else:
        time_pos = len(model.t) - 1
    return time_pos


def check_dictionary(value):
    """Check if variable is a dictionary, and raise an error if it is not."""
    if not isinstance(value, dict):
        raise TypeError(f"Input must be a dict, but is {type(value)} instead.")


class DomainValueError(Exception):
    """Exception raised for values that are outside their possible domain.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, parameter, expected_domain, unexpected_value, prefix=""):
        """Initialize error class.

        Args:
            parameter (str): Name of invalid parameter.
            expected_domain (str): Description of domain that was expected of parameter.
            unexpected_value (numerical): Value of parameter that is not in expected domain.
            prefix (str): Additional error description to be put before standardized error message.
        """
        if len(prefix) and prefix[-1] != " ":
            prefix += " "
        self.message = f"{prefix}{parameter} must be {expected_domain}, but was {unexpected_value} instead."
        super().__init__(self.message)


class MissingValueError(Exception):
    """Exception raised when one or more required parameters are missing.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        """Initialize error class."""
        self.message = message
        super().__init__(self.message)
