"""Author: Jorrit Bakker.

File testing functionality of mass_balance module.
"""

import pytest
from mibitrans.analysis.mass_balance import mass_balance
from mibitrans.data.parameters import ModelParameters
from tests.test_example_data import testing_massbalance_instant
from tests.test_example_data import testing_massbalance_lindecay
from tests.test_example_data import testing_massbalance_nodecay


@pytest.fixture(scope="module")
def test_model_pars():
    """ModelParameters fixture with increased spatial resolution, specifically for testing mass balance."""
    return ModelParameters(model_length=50, model_width=30, model_time=3 * 365, dx=1, dy=1, dt=1 * 365)


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_domenico_nodecay_model", testing_massbalance_nodecay),
        ("test_domenico_lineardecay_model", testing_massbalance_lindecay),
        ("test_domenico_instantreaction_model", testing_massbalance_instant),
        (ModelParameters(), TypeError),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_balance_numerical(model, expected, request) -> None:
    """Test if mass balance is correctly calculated by comparing to precomputed results."""
    if isinstance(expected, dict):
        model_object = request.getfixturevalue(model)
        dictionary = mass_balance(model_object, time=3 * 365)
        for key, output_item in dictionary.items():
            assert expected[key] == pytest.approx(output_item)
    else:
        with pytest.raises(expected):
            mass_balance(model, time=3 * 365)
