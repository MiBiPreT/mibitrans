"""Author: Jorrit Bakker.

Module handling data input in the form of a dictionary.
"""
import warnings
import numpy as np
from mibitrans.data.parameter_information import key_dictionary, electron_acceptor_utilization
from mibitrans.data.check_input import (_check_float_positive, _check_float_fraction, _check_float_retardation,
                                        _check_array_float_positive, _check_total_mass)
from dataclasses import dataclass


@dataclass
class HydrologicalParameters:
    """Dataclass handling input of hydrological parameters.

    Args:
        velocity (float) : Flow velocity in the direction of the groundwater gradient, in [m/d]. Optional if h_gradient and h_conductivity are specified.
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
    velocity : float = None
    h_gradient : float = None
    h_conductivity : float = None
    porosity : float = None
    alpha_x : float = None
    alpha_y : float = None
    alpha_z : float = 1e-10
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

        if self.verbose:
            print("All required hydrological input arguments are present.")

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

        if self.verbose:
            print("All hydrological input arguments are valid.")

        # Velocity is calculated from hydraulic gradient and conductivity when both are given.
        if self.h_gradient and self.h_conductivity:
            # Giving h_gradient and h_conductivity is more specific than giving velocity. So input velocity will be overridden.
            if self.velocity is not None:
                warnings.warn("Both velocity and h_gradient & h_conductivity are defined. Value for velocity will be overridden.", UserWarning)
            self.velocity = self.h_gradient * self.h_conductivity / self.porosity
            if self.verbose:
                print(f"Groundwater flow velocity has been calculated to be {self.velocity} m/d.")

@dataclass
class AdsorptionParameters:
    """Dataclass handling adsorption parameters.

    Args:
        retardation (float) : Retardation factor for transported contaminant [-]. Optional if bulk_density,
            partition_coefficient and fraction_organic_carbon are specified.
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
    retardation : float = None
    bulk_density : float = None
    partition_coefficient : float = None
    fraction_organic_carbon : float = None
    verbose : bool = False

    def __post_init__(self):
        if self.retardation is None and (self.bulk_density is None
                                         or self.bulk_density is None
                                         or self.partition_coefficient is None
                                         or self.fraction_organic_carbon is None):
            raise ValueError("AdsorptionParameters missing required arguments: either retardation or (bulk_density, partition_coefficient and fraction_organic_carbon).")

        if self.verbose:
            print("All required adsorption input arguments are present.")

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

        if self.verbose:
            print("All adsorption input arguments are valid.")

    # Retardation factor is not calculated from bulk density, partition coefficient and fraction organic carbon
    # in this dataclass, since porosity is required as well, which is defined in the HydrologicalParameters class.

    def calculate_retardation(self, porosity : float):
        """Calculate retardation factor from other input if not given"""
        if self.retardation is None:
            self.retardation = 1 + (self.bulk_density / porosity) * self.partition_coefficient * self.fraction_organic_carbon
            if self.verbose:
                print(f"Retardation factor has been calculated to be {self.retardation}.")

# Model extent in the x direction (length) in [m]
# Model extent in the y direction (width) in [m]
# Model end time in [years]
# Thickness (z-direction, depth) of source zone in [m]
# Concentration in the source zone as array with [y-location, concentration] in [[m], [g/m^3]]
# Mass of contaminant source in [kg], or "inf" for infinite source.

@dataclass
class DegradationParameters:
    """Dataclass handling degradation parameters.

    Args:
        decay_rate (float) : First order (linear) decay coefficient in [1/day]. Only required for linear decay models.
            Optional if half_life is specified.
        half_life (float) : Contaminant half life for 1st order (linear) decay, in [days]. Only required for
            linear decay models. Optional if half_life is specified.
        delta_oxygen (float) : Difference between background oxygen and plume oxygen concentrations, in [g/m^3].
            Only required for instant reaction models.
        delta_nitrate (float) : Difference between background nitrate and contaminant plume nitrate concentrations,
            in [g/m^3]. Only required for instant reaction models.
        ferrous_iron (float) : Ferrous iron concentration in contaminant plume, in [g/m^3]. Only required for
            instant reaction models.
        delta_sulfate (float) : Difference between background sulfate and plume sulfate concentrations, in [g/m^3].
            Only required for instant reaction models.
        methane (float) : Methane concentration in contaminant plume, in [g/m^3]. Only required for
            instant reaction models.
        verbose (bool, optional): Verbose mode. Defaults to False.

    Methods:
        utilization_factor : Customize electron acceptor utilization factors. By default, electron acceptor utilization
            factors for a BTEX mixture are used, based on values by Wiedemeier et al. (1995), see
            electron_acceptor_utilization in mibitrans.data.parameter_information.

    Raises:
        ValueError : If input parameters are incomplete or outside the valid domain.
        TypeError : If input parameters of incorrect datatype.

    """
    decay_rate : float = None
    half_life : float = None
    delta_oxygen : float = None
    delta_nitrate : float = None
    ferrous_iron : float = None
    delta_sulfate : float = None
    methane : float = None
    verbose : bool = False

    def __post_init__(self):
        if (self.decay_rate is None and self.half_life is None) and (self.delta_oxygen is None
                                                                     or self.delta_nitrate is None
                                                                     or self.ferrous_iron is None
                                                                     or self.delta_sulfate is None
                                                                     or self.methane is None):
            raise ValueError("DegradationParameters missing missing required arguments: either decay rate or half life, or electron acceptor/donor concentrations.")

        # All input parameters should be positive floats, raise appropriate error if not
        for parameter, value in self.__dict__.items():
            if parameter == "verbose":
                continue
            else:
                error = _check_float_positive(parameter, value)

            if error and (value is not None):
                raise error

        if self.half_life:
            decay_rate = np.log(2) / self.half_life
            if self.decay_rate and (self.decay_rate != decay_rate):
                warnings.warn("Both contaminant decay rate constant and half life are defined, but are not equal. Only value for decay rate constant will be used in calculations.", UserWarning)
            else:
                self.decay_rate = decay_rate

        self.electron_acceptor_utilization = electron_acceptor_utilization

    def require_electron_acceptor(self):
        missing_ea = []
        if self.delta_oxygen is None:
            missing_ea.append("delta_oxygen")
        if self.delta_nitrate is None:
            missing_ea.append("delta_nitrate")
        if self.ferrous_iron is None:
            missing_ea.append("ferrous_iron")
        if self.delta_sulfate is None:
            missing_ea.append("delta_sulfate")
        if self.methane is None:
            missing_ea.append("methane")

        if len(missing_ea) > 0:
            raise ValueError(f"Instant reaction model requires concentrations of {missing_ea}.")

    def require_linear_decay(self):
        if self.decay_rate is None and self.half_life is None:
            raise ValueError("Linear reaction model requires decay rate or half life.")

    def utilization_factor(self,
                           util_oxygen : float = None,
                           util_nitrate : float = None,
                           util_ferrous_iron : float = None,
                           util_sulfate : float = None,
                           util_methane : float = None
                           ):
        """Introduce custom utilization factors for each electron donor/acceptor species. By default, utilization factors for mix of BTEX are used.

        Args:
            util_oxygen (float) : utilization factor of oxygen, as mass of oxygen consumed per mass of biodegradated contaminant [g/g].
            util_nitrate (float) : utilization factor of nitrate, as mass of nitrate consumed per mass of biodegrated contaminant [g/g].
            util_ferrous_iron (float) : utilization factor of ferrous iron, as mass of ferrous iron generated per mass of biodegradated contaminant [g/g].
            util_sulfate (float) : utilization factor of sulfate, as mass of sulfate consumed per mass of biodegradated contaminant [g/g].
            util_methane (float) : utilization factor of methane, as mass of methane generated per mass of biodegrated contaminant [g/g].

        Raises:
            ValueError : If input parameters are incomplete or outside the valid domain.
            TypeError : If input parameters of incorrect datatype.

        """
        utils = locals()
        for parameter, value in utils.items():
            if parameter != "self" and value is not None:
                error = _check_float_positive(parameter, value)
                if error:
                    raise error
                self.electron_acceptor_utilization[parameter] = value

@dataclass
class SourceParameters:
    """Dataclass handling source parameters. Specifying concentrations and extent of source zone.

    Args:
        source_zone_boundary (np.ndarray) : Outer boundary of each source zone, in transverse horizontal direction (y-coordiante) [m].
            y=0 is at the middle of the contaminant source. Input as numpy array of length equal to the amount of source zone.
            Last value in the array is the limit of the source. For a source with a single source zone, only one value is required.
            Source is symmetrical in the x-axis.
        source_zone_concentration (np.ndarray) : Contaminant concentration in each source zone [g/m^3]. Input as numpy
            array in the same order and of the same length as specified in source_zone_boundary.
        depth (float) : Depth (transeverse vertical or z-dimension) of the source zone in [m].
        total_mass (float | str) : Mass of contaminant present in source zone, either expressed in [g], or set to 'infinite'.
            The latter meaning that the source mass and therefore, the source zone concentrations do not diminish over time.
        verbose (bool, optional): Verbose mode. Defaults to False.

    Raises:
        ValueError : If input parameters are incomplete or outside the valid domain.
        TypeError : If input parameters of incorrect datatype.
    """

    source_zone_boundary : np.ndarray = None
    source_zone_concentration : np.ndarray = None
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
    """Dataclass handling model discretization parameters.

    Args:
        model_length (float) : Model extent in the longitudinal (x) direction in [m].
        model_width (float) : Model extent in the transverse horizontal (y) direction in [m].
        model_time (float) : Model duration in [days].
        dx (float) : Model grid discretization step size in the longitudinal (x) direction, in [m].
        dy (float) : Model grid discretization step size in the transverse horizontal (y) direction, in [m].
        dt (float) : Model time discretization step size, in [days]
        verbose (bool, optional): Verbose mode. Defaults to False.

    Raises:
        ValueError : If input parameters are incomplete or outside the valid domain.
        TypeError : If input parameters of incorrect datatype.

    """
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
