"""Author: Jorrit Bakker.

File testing functionality of mass_balance module.
"""

import pytest
import anatrans.analysis.mass_balance as mb
from anatrans.data.parameter_information import testing_dictionary
from tests.test_example_results import testing_massbalance_instant
from tests.test_example_results import testing_massbalance_lindecay
from tests.test_example_results import testing_massbalance_nodecay


@pytest.mark.parametrize(
    "input, time, dt, expected",
    [
        (testing_dictionary, 2, 1, 2),
        (testing_dictionary, None, None, 3),
        (testing_dictionary, 6, 1, 3),
        (testing_dictionary, 1.7, 1, 2),
    ])

def test_balance_time(input : dict, time, dt, expected) -> None:
    """Tests time point determination of balance function."""
    obj_mb = mb.MassBalance(input, dx=None, dy=None, dt=dt, mode="no_decay")
    output = obj_mb.balance(time=time)
    assert output["time"] == pytest.approx(expected)

@pytest.mark.parametrize(
    "input, stepsize, time, mode, expected",
    [
        (testing_dictionary, (1,1,1), 3, "no_decay", testing_massbalance_nodecay),
        (testing_dictionary, (1,1,1), 3, "linear_decay", testing_massbalance_lindecay),
        (testing_dictionary, (1,1,1), 3, "instant_reaction", testing_massbalance_instant),
    ])

def test_balance_results(input : dict, stepsize, time, mode : str, expected) -> None:
    """Test if mass balance is correctly calculated by comparing to precomputed results."""
    dx, dy, dt = stepsize
    obj_mb = mb.MassBalance(input, dx=dx, dy=dy, dt=dt, mode=mode)
    output = obj_mb.balance(time=time)
    for key, output_item in output.items():
        assert expected[key] == pytest.approx(output_item)
