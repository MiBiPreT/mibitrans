"""Author: Jorrit Bakker.

File containing various dictionaries used for evaluation of names, value types and units of input data.

"""

from dataclasses import dataclass
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
    """Make UtilizationFactor object.

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

    def __post_init__(self):
        """Initialize utilization factors as dictionary."""
        self.dictionary = dict(
            util_oxygen=self.util_oxygen,
            util_nitrate=self.util_nitrate,
            util_ferrous_iron=self.util_ferrous_iron,
            util_sulfate=self.util_sulfate,
            util_methane=self.util_methane,
        )
