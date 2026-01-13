"""Author: Jorrit Bakker.

File containing various dictionaries used for evaluation of names, value types and units of input data.

"""

from dataclasses import dataclass
import numpy as np
from mibitrans.data.check_input import validate_input_values

# Couples utilization factors to electron acceptors/donors
util_to_conc_name = {
    "util_oxygen": "delta_oxygen",
    "util_nitrate": "delta_nitrate",
    "util_ferrous_iron": "ferrous_iron",
    "util_sulfate": "delta_sulfate",
    "util_methane": "methane",
}

mass_balance_renaming_dictionary = {
    "source_mass_0": "mass t = 0",
    "source_mass_t": "mass t = ",
    "source_mass_change": "delta mass",
    "plume_mass_no_decay": "plume mass",
    "transport_outside_extent_nodecay": "mass transported outside model extent",
    "plume_mass_linear_decay": "plume mass",
    "transport_outside_extent_lineardecay": "mass transported outside model extent",
    "plume_mass_degraded_linear": "plume mass degraded",
    "source_mass_instant_t": "source mass t = ",
    "source_mass_instant_change": "delta source mass",
    "plume_mass_no_decay_instant_reaction": "plume mass before decay",
    "plume_mass_instant_reaction": "plume mass after decay",
    "plume_mass_degraded_instant": "plume mass degraded",
    "electron_acceptor_mass_change": "change in mass (kg)",
}


@dataclass
class UtilizationFactor:
    """Make object containing information about electron acceptor utilization factor.

    Args:
        util_oxygen (float) : utilization factor of oxygen, as mass of oxygen consumed
            per mass of biodegraded contaminant [g/g].
        util_nitrate (float) : utilization factor of nitrate, as mass of nitrate consumed
            per mass of biodegraded contaminant [g/g].
        util_ferrous_iron (float) : utilization factor of ferrous iron, as mass of ferrous iron generated
            per mass of biodegraded contaminant [g/g].
        util_sulfate (float) : utilization factor of sulfate, as mass of sulfate consumed
            per mass of biodegraded contaminant [g/g].
        util_methane (float) : utilization factor of methane, as mass of methane generated
            per mass of biodegraded contaminant [g/g].

    Raises:
        ValueError : If input parameters are incomplete or outside the valid domain.
        TypeError : If input parameters of incorrect datatype.

    """

    util_oxygen: float
    util_nitrate: float
    util_ferrous_iron: float
    util_sulfate: float
    util_methane: float

    def __setattr__(self, parameter, value):
        """Override parent method to validate input when attribute is set."""
        if parameter != "dictionary":
            validate_input_values(parameter, value)
        super().__setattr__(parameter, value)

    @property
    def dictionary(self):
        """Returns utilization factors in the form of a dictionary."""
        return dict(
            util_oxygen=self.util_oxygen,
            util_nitrate=self.util_nitrate,
            util_ferrous_iron=self.util_ferrous_iron,
            util_sulfate=self.util_sulfate,
            util_methane=self.util_methane,
        )

    @property
    def array(self):
        """Return utilization factors in the form of an array, in order of [O2, NO3, Fe, SO4, CH4]."""
        return np.array(
            [self.util_oxygen, self.util_nitrate, self.util_ferrous_iron, self.util_sulfate, self.util_methane]
        )


@dataclass
class ElectronAcceptors:
    """Make object with concentrations of electron acceptors.

    Dataclass which handles the entry of electron acceptor concentrations used for the instant reaction biodegradation
    method. As plume concentrations for reduced electron acceptor species and as difference between plume and background
    concentrations for the electron acceptors themselves.

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
    """

    delta_oxygen: float
    delta_nitrate: float
    ferrous_iron: float
    delta_sulfate: float
    methane: float

    def __setattr__(self, parameter, value):
        """Override parent method to validate input when attribute is set."""
        validate_input_values(parameter, value)
        super().__setattr__(parameter, value)

    @property
    def dictionary(self):
        """Returns electron acceptors in the form of a dictionary."""
        return dict(
            delta_oxygen=self.delta_oxygen,
            delta_nitrate=self.delta_nitrate,
            ferrous_iron=self.ferrous_iron,
            delta_sulfate=self.delta_sulfate,
            methane=self.methane,
        )

    @property
    def array(self):
        """Return electron acceptor concentrations in the form of an array, in order of [O2, NO3, Fe, SO4, CH4]."""
        return np.array([self.delta_oxygen, self.delta_nitrate, self.ferrous_iron, self.delta_sulfate, self.methane])
