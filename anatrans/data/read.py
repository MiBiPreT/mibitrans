"""Author: Jorrit Bakker.

Module handling data input in the form of a dictionary.
"""
import sys

sys.path.insert(1, "./anatrans/data")

from parameter_information import key_dictionary


def from_dict(dictionary: dict,
              check_input: bool = True,
              verbose: bool = True
              ) -> dict:
    """Format and structure input dictionary into a standardized dictionary."""
    params = {}
    unknown_keys = []
    for key_input, value_input in dictionary.items():
        key_in_known_keys = False
        for key_params, key_known in key_dictionary.items():
            if key_input in key_known:
                params[key_params] = value_input
                key_in_known_keys = True
                break
        if not key_in_known_keys:
            unknown_keys.append(key_input)

    if verbose and len(unknown_keys) > 0:
        print("The following keys were not recognized and not included in output dictionary:", unknown_keys)
    elif verbose:
        print("All keys were recognized")
    return params
