"""Author: Jorrit Bakker.

File testing functionality of mass_balance module.
"""

import os
import numpy as np
import pytest
from mibitrans.analysis.mass_balance import mass_balance
from mibitrans.data.parameters import ModelParameters
from mibitrans.data.parameters import SourceParameters
from mibitrans.transport.models import Anatrans
from mibitrans.transport.models import Mibitrans
from tests.test_example_data import testing_massbalance_instant_dom
from tests.test_example_data import testing_massbalance_instant_dom_inf
from tests.test_example_data import testing_massbalance_instant_kar
from tests.test_example_data import testing_massbalance_instant_kar_inf
from tests.test_example_data import testing_massbalance_lindecay_dom
from tests.test_example_data import testing_massbalance_lindecay_kar
from tests.test_example_data import testing_massbalance_nodecay_dom
from tests.test_example_data import testing_massbalance_nodecay_kar

IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"


@pytest.fixture(scope="module")
def test_model_pars():
    """ModelParameters fixture with increased spatial resolution, specifically for testing mass balance."""
    return ModelParameters(model_length=50, model_width=30, model_time=3 * 365, dx=1, dy=1, dt=1 * 365)


@pytest.fixture(scope="module")
def test_source_pars_inf():
    """SourceParameters fixture with example data for tests."""
    return SourceParameters(
        source_zone_boundary=np.array([5, 10, 15]),
        source_zone_concentration=np.array([10, 5, 2]),
        depth=10,
        total_mass="inf",
    )


@pytest.fixture(scope="module")
def test_mibitrans_nodecay_model_mb(test_hydro_pars, test_att_pars_nodecay, test_source_pars, test_model_pars):
    """Mibitrans with no decay fixture model object for testing."""
    obj = Mibitrans(test_hydro_pars, test_att_pars_nodecay, test_source_pars, test_model_pars)
    obj.run()
    return obj


@pytest.fixture(scope="module")
def test_mibitrans_lineardecay_model_mb(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """Mibitrans with linear decay fixture model object for testing."""
    obj = Mibitrans(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    obj.run()
    return obj


@pytest.fixture(scope="module")
def test_mibitrans_instantreaction_model_mb(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """Mibitrans with instant reacion fixture model object for testing."""
    obj = Mibitrans(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    obj.instant_reaction(dict(delta_oxygen=0.5, delta_nitrate=0.5, ferrous_iron=0.5, delta_sulfate=0.5, methane=0.5))
    obj.run()
    return obj


@pytest.fixture(scope="module")
def test_mibitrans_instantreaction_model_mb_inf(test_hydro_pars, test_att_pars, test_source_pars_inf, test_model_pars):
    """Mibitrans with instant reaction and infinite source mass fixture model object for testing."""
    obj = Mibitrans(test_hydro_pars, test_att_pars, test_source_pars_inf, test_model_pars)
    obj.instant_reaction(dict(delta_oxygen=0.5, delta_nitrate=0.5, ferrous_iron=0.5, delta_sulfate=0.5, methane=0.5))
    obj.run()
    return obj


@pytest.fixture(scope="module")
def test_anatrans_instantreaction_model_inf(test_hydro_pars, test_att_pars, test_source_pars_inf, test_model_pars):
    """Anatrans instant reaction fixture model object for testing."""
    obj = Anatrans(test_hydro_pars, test_att_pars, test_source_pars_inf, test_model_pars)
    obj.instant_reaction(dict(delta_oxygen=0.5, delta_nitrate=0.5, ferrous_iron=0.5, delta_sulfate=0.5, methane=0.5))
    obj.run()
    return obj


@pytest.mark.skipif(IN_GITHUB_ACTIONS, reason="Mass balance tests temporarily disabled due to incoming refactor.")
@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_mibitrans_nodecay_model_mb", testing_massbalance_nodecay_kar),
        ("test_mibitrans_lineardecay_model_mb", testing_massbalance_lindecay_kar),
        ("test_mibitrans_instantreaction_model_mb", testing_massbalance_instant_kar),
        ("test_mibitrans_instantreaction_model_mb_inf", testing_massbalance_instant_kar_inf),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_balance_numerical_mibitrans(model, expected, request) -> None:
    """Test if mass balance is correctly calculated by comparing to precomputed results for Mibitrans model."""
    model_object = request.getfixturevalue(model)
    dictionary = mass_balance(model_object, time=3 * 365)
    for key, output_item in dictionary.items():
        assert expected[key] == pytest.approx(output_item)


@pytest.mark.skipif(IN_GITHUB_ACTIONS, reason="Mass balance tests temporarily disabled due to incoming refactor.")
@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_anatrans_model_nodecay", testing_massbalance_nodecay_dom),
        ("test_anatrans_model_lineardecay", testing_massbalance_lindecay_dom),
        ("test_anatrans_model_instantreaction", testing_massbalance_instant_dom),
        ("test_anatrans_instantreaction_model_inf", testing_massbalance_instant_dom_inf),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_balance_numerical(model, expected, request) -> None:
    """Test if mass balance is correctly calculated by comparing to precomputed results for Anatrans model."""
    model_object = request.getfixturevalue(model)
    dictionary = mass_balance(model_object, time=3 * 365)
    for key, output_item in dictionary.items():
        assert expected[key] == pytest.approx(output_item)
