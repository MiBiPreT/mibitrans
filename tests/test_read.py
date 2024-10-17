"""Author: Jorrit Bakker.

Module handling data input in the form of a dictionary.
"""

import pytest
import anatrans.data.read as rd
from anatrans.data.parameter_names import fromdict_dictionary as f_dict

@pytest.mark.parametrize(
    "test, expected",
    [
        ({f_dict["v"][0] : 1, f_dict["R"][0] : 2, f_dict["mu"][0] : 3}, dict(v = 1, R = 2, mu = 3)),
        (dict(v = 1, R = 2, nonsense = 3), dict(v = 1, R = 2)),
        (dict(nonsense = 3), dict()),
        (dict(), dict()),
        ({f_dict["v"][-1] : 1}, {"v" : 1})
    ])

def test_from_dict(test, expected) -> None:
    """Test if from_dict gives expected output for various input dictionaries."""
    result = rd.from_dict(test)
    assert result == expected
