"""Author: Jorrit Bakker.

Module evaluating if a dictionary contains all required (correct) parameters for analysis
"""

import warnings
import numpy as np
import mibitrans


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
            return DomainValueError(f"{parameter} must be >= 0")
    else:
        return is_float


def _check_numeric_fraction(parameter: str, value):
    """Check if a variable is numerical and if it is between 0 and 1."""
    is_float = _check_numeric(parameter, value)
    if is_float is None:
        if 0 <= value <= 1:
            return None
        else:
            return DomainValueError(f"{parameter} must be between 0 and 1")
    else:
        return is_float


def _check_numeric_retardation(parameter: str, value):
    """Check if a variable is numerical and if it is 1 or larger."""
    is_float = _check_numeric(parameter, value)
    if is_float is None:
        if value >= 1:
            return None
        else:
            return DomainValueError(f"{parameter} must be 1 or larger.")
    else:
        return TypeError(f"{parameter} must be a float, but is {type(value)} instead.")


def _check_list_array_positive(parameter: str, value, sublist_allowed):
    """Check if variable contains numeric values or lists/arrays of numeric values if allowed."""
    if all(isinstance(item, (list, np.ndarray)) for item in value) and sublist_allowed:
        if isinstance(value, np.ndarray):
            if len(value.shape) > 1:
                return ValueError(
                    f"{parameter} should either be a 1-dimensional array, a list of 1-dimensional arrays "
                    f"or list of lists, but is {len(value.shape)}-dimensional array instead."
                )
        for item in value:
            if not all(isinstance(elem, (int, float, np.floating, np.integer)) for elem in item):
                return TypeError(f"All sub-elements of {parameter} should be a float.")
            if not all(elem >= 0 for elem in item):
                return DomainValueError(f"All sub-elements in {parameter} should be >= 0.")
        return None
    elif all(isinstance(item, (int, float, np.floating, np.integer)) for item in value):
        if all(item >= 0 for item in value):
            return None
        else:
            return DomainValueError(f"All elements of {parameter} should be >= 0.")
    else:
        if isinstance(value, np.ndarray):
            if len(value.shape) != 1:
                return ValueError(f"{parameter} must be a 1D array of floats or list of 1D-arrays/floats.")
        if sublist_allowed:
            return TypeError(f"All elements of {parameter} should either be a float or list/array of floats.")
        else:
            return TypeError(f"All elements of {parameter} should be a float.")


def _check_array_list_numeric_positive(parameter: str, value, sublist_allowed):
    """Check if variable is numpy array, list, or numerical, if it is positive and if an array is 1-dimensional."""
    if isinstance(value, (np.ndarray, list)):
        return _check_list_array_positive(parameter, value, sublist_allowed)
    elif isinstance(value, (float, int, np.floating, np.integer)):
        return _check_numeric_positive(parameter, value)
    else:
        return TypeError(f"{parameter} must be a float, list or numpy array, but is {type(value)} instead.")


def validate_source_zones(boundary, concentration):
    """Validate and adapt input of source_zone_boundary and source_zone_concentration arrays."""
    chain_exception = False
    # Ensure boundary and concentration are numpy arrays
    if isinstance(boundary, (float, int, np.floating, np.integer)):
        boundary = np.array([boundary])
    else:
        boundary = np.array(boundary)

    if isinstance(concentration, list):
        if isinstance(concentration[0], (list, np.ndarray)):
            if len(concentration) == 1:
                concentration = concentration[0]
            else:
                concentration = [np.array(conc) for conc in concentration]
        else:
            if len(boundary) != len(concentration) and len(boundary) == 1:
                # When only a single source zone boundary is given, but multiple (single) source zone concentrations,
                # it is interpreted as varying single source concentrations for purpose of chain decay. Therefore, no
                # error will be raised if length of boundary array != length concentration array.
                chain_exception = True
                concentration = [np.array([conc]) for conc in concentration]
            else:
                concentration = np.array(concentration)
    elif isinstance(concentration, np.ndarray):
        if len(boundary) != len(concentration) and len(boundary) == 1:
            chain_exception = True
    else:
        concentration = np.array([concentration])

    # Reording of source zone boundary if not in correct order decrepit from v1.1.0 due to conflict with chain decay
    # Furthermore, it is considered good practice to purposefully set source zone boundary in correct order. To prevent
    # unintended source discretization.
    if len(boundary) > 1:
        if not all(boundary[:-1] <= boundary[1:]):
            boundary.sort()
            raise ValueError(
                "source_zone_boundary locations should be ordered by distance from source zone center. Thus, current "
                f"input for source_zone_boundary is supposed to be {boundary}. source_zone_concentration should be re-"
                f"ordered accordingly as well, with highest concentrations at the innermost source zone."
            )
    # Superposition method only works if the zone closer to the center has higher concentration than outer zones
    check_conc = [concentration] if not isinstance(concentration, list) else concentration
    for conc in check_conc:
        # Each given source zone boundary should have a corresponding concentration, and vice versa
        if boundary.shape != conc.shape and not chain_exception:
            try:
                len_conc = len(conc)
            except TypeError:
                conc = np.array([conc])
                len_conc = len(conc)

            raise ValueError(
                f"Length of source zone boundary (len={len(boundary)}, for {boundary}) and source zone concentration "
                f"(len={len_conc}, for {conc}) do not match. Make sure they are of equal length."
            )
        if (
            not all(conc[:-1] >= conc[1:])
            and not isinstance(conc, (float, int, np.floating, np.integer))
            and not chain_exception
        ):
            raise ValueError(
                "Source zone concentrations should be in descending order; no source zone can have a concentration "
                "higher than the concentration of a zone closer to source center, due to the superposition method."
            )
    return boundary, concentration


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
            return DomainValueError(f"{parameter} must be >= 0, or set to 'infinite'.")
    else:
        return TypeError(f"{parameter} must be a float or 'infinite', but is {type(value)} instead.")


def check_dictionary(value):
    """Check if variable is a dictionary, and raise an error if it is not."""
    if not isinstance(value, dict):
        raise TypeError(f"Input must be a dict, but is {type(value)} instead.")


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


def validate_input_values(parameter, value):
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
        case "electron_acceptor":
            error = _check_dataclass(parameter, value, mibitrans.data.parameter_information.FringeElectronAcceptors)
        # Parameters which can be any float value
        case "y_position":
            error = _check_numeric(parameter, value)
        # Parameters which have domain [0,1]
        case "porosity" | "fraction_organic_carbon":
            error = _check_numeric_fraction(parameter, value)
        # Parameters which are input as single values, lists or numpy arrays
        case (
            "source_zone_boundary"
            | "decay_rate"
            | "half_life"
            | "mass_ratios"
            | "electron_acceptor_concentration"
            | "stoichiometric_ratio"
            | "molecular_weight_electron_acceptor"
        ):
            error = _check_array_list_numeric_positive(parameter, value, sublist_allowed=False)
        case "source_zone_concentration":
            error = _check_array_list_numeric_positive(parameter, value, sublist_allowed=True)
        case "electron_acceptors":
            error = _check_electron_acceptor(value)
        # All other parameters are checked as floats on positive domain
        case _:
            error = _check_numeric_positive(parameter, value)

    if error and (value is not None):
        raise error


class DomainValueError(Exception):
    """Exception raised for values that are outside their possible domain.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        """Initialize error class."""
        self.message = message
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
