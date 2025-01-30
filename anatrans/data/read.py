"""Author: Jorrit Bakker.

Module handling data input in the form of a dictionary.
"""
from anatrans.data.parameter_information import key_dictionary


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
