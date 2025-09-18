"""Author: Jorrit Bakker.

File testing functionality of mass_balance module.
"""

import pytest
from mibitrans.analysis.mass_balance import mass_balance
from mibitrans.transport.domenico import InstantReaction
from mibitrans.transport.domenico import LinearDecay
from mibitrans.transport.domenico import NoDecay
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
        (NoDecay(test_hydro_pars, test_ads_pars, test_source_pars, test_model_pars), testing_massbalance_nodecay),
        (
            LinearDecay(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars),
            testing_massbalance_lindecay,
        ),
        (
            InstantReaction(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars),
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
