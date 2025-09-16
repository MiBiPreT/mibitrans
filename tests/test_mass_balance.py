"""Author: Jorrit Bakker.

File testing functionality of mass_balance module.
"""

import pytest
import mibitrans.analysis.mass_balance as mb
from mibitrans.analysis.mass_balance import mass_balance
from mibitrans.data.parameter_information import testing_dictionary
from mibitrans.transport.domenico import instant_reaction
from mibitrans.transport.domenico import linear_decay
from mibitrans.transport.domenico import no_decay
from tests.test_example_data import test_ads_pars
from tests.test_example_data import test_deg_pars
from tests.test_example_data import test_hydro_pars
from tests.test_example_data import test_model_pars
from tests.test_example_data import test_source_pars
from tests.test_example_data import testing_massbalance_instant
from tests.test_example_data import testing_massbalance_lindecay
from tests.test_example_data import testing_massbalance_nodecay

test_model_pars.dx = 1
test_model_pars.dy = 1
test_model_pars.dt = 1


@pytest.mark.parametrize(
    "model, expected",
    [
        (no_decay(test_hydro_pars, test_ads_pars, test_source_pars, test_model_pars), testing_massbalance_nodecay),
        (
            linear_decay(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars),
            testing_massbalance_lindecay,
        ),
        (
            instant_reaction(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars),
            testing_massbalance_instant,
        ),
        (test_hydro_pars, TypeError),
    ],
)
def test_balance_numerical(model, expected) -> None:
    """Test if mass balance is correctly calculated by comparing to precomputed results."""
    if isinstance(expected, dict):
        dictionary = mass_balance(model, time=3 * 365)
        for key, output_item in dictionary.items():
            assert expected[key] == pytest.approx(output_item)
    else:
        with pytest.raises(expected):
            mass_balance(model, time=3 * 365)


########################################################################################################################
####################################### Pre-refactor functionalities, decrepit #########################################
########################################################################################################################


@pytest.mark.parametrize(
    "input, time, dt, expected",
    [
        (testing_dictionary, 2, 1, 2),
        (testing_dictionary, None, None, 3),
        (testing_dictionary, 6, 1, 3),
        (testing_dictionary, 1.7, 1, 2),
    ],
)
def test_balance_time(input: dict, time, dt, expected) -> None:
    """Tests time point determination of balance function."""
    obj_mb = mb.MassBalance(input, dx=None, dy=None, dt=dt, mode="no_decay")
    output = obj_mb.balance(time=time)
    assert output["time"] == pytest.approx(expected)


@pytest.mark.parametrize(
    "input, stepsize, time, mode, expected",
    [
        (testing_dictionary, (1, 1, 1), 3, "no_decay", testing_massbalance_nodecay),
        (testing_dictionary, (1, 1, 1), 3, "linear_decay", testing_massbalance_lindecay),
        (testing_dictionary, (1, 1, 1), 3, "instant_reaction", testing_massbalance_instant),
    ],
)
def test_balance_results(input: dict, stepsize, time, mode: str, expected) -> None:
    """Test if mass balance is correctly calculated by comparing to precomputed results."""
    expected["time"] = expected["time"] / 365
    dx, dy, dt = stepsize
    obj_mb = mb.MassBalance(input, dx=dx, dy=dy, dt=dt, mode=mode)
    output = obj_mb.balance(time=time)
    for key, output_item in output.items():
        assert expected[key] == pytest.approx(output_item)
