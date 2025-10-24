"""Author: Jorrit Bakker.

File testing functionality of mass_balance module.
"""

import numpy as np
import pytest
from mibitrans.analysis.mass_balance import mass_balance
from mibitrans.data.parameters import ModelParameters
from mibitrans.data.parameters import SourceParameters
from mibitrans.transport import domenico as dom
from mibitrans.transport import karanovic as kar
from tests.test_example_data import testing_massbalance_instant_dom
from tests.test_example_data import testing_massbalance_instant_dom_inf
from tests.test_example_data import testing_massbalance_instant_kar
from tests.test_example_data import testing_massbalance_instant_kar_inf
from tests.test_example_data import testing_massbalance_lindecay_dom
from tests.test_example_data import testing_massbalance_lindecay_kar
from tests.test_example_data import testing_massbalance_nodecay_dom
from tests.test_example_data import testing_massbalance_nodecay_kar


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
def test_karanovic_nodecay_model_mb(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """karanovic.NoDecay fixture model object for testing."""
    return kar.NoDecay(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)


@pytest.fixture(scope="module")
def test_karanovic_lineardecay_model_mb(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """karanovic.LinearDecay fixture model object for testing."""
    return kar.LinearDecay(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)


@pytest.fixture(scope="module")
def test_karanovic_instantreaction_model_mb(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """karanovic.InstantReaction fixture model object for testing."""
    return kar.InstantReaction(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)


@pytest.fixture(scope="module")
def test_karanovic_instantreaction_model_mb_inf(test_hydro_pars, test_att_pars, test_source_pars_inf, test_model_pars):
    """karanovic.InstantReaction fixture model object for testing."""
    return kar.InstantReaction(test_hydro_pars, test_att_pars, test_source_pars_inf, test_model_pars)


@pytest.fixture(scope="module")
def test_domenico_instantreaction_model_inf(test_hydro_pars, test_att_pars, test_source_pars_inf, test_model_pars):
    """karanovic.InstantReaction fixture model object for testing."""
    return dom.InstantReaction(test_hydro_pars, test_att_pars, test_source_pars_inf, test_model_pars)


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_karanovic_nodecay_model_mb", testing_massbalance_nodecay_kar),
        ("test_karanovic_lineardecay_model_mb", testing_massbalance_lindecay_kar),
        ("test_karanovic_instantreaction_model_mb", testing_massbalance_instant_kar),
        ("test_karanovic_instantreaction_model_mb_inf", testing_massbalance_instant_kar_inf),
        (ModelParameters(), TypeError),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_balance_numerical_karanovic(model, expected, request) -> None:
    """Test if mass balance is correctly calculated by comparing to precomputed results for Karanovic model."""
    if isinstance(expected, dict):
        model_object = request.getfixturevalue(model)
        dictionary = mass_balance(model_object, time=3 * 365)
        for key, output_item in dictionary.items():
            assert expected[key] == pytest.approx(output_item)
    else:
        with pytest.raises(expected):
            mass_balance(model, time=3 * 365)


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_domenico_nodecay_model", testing_massbalance_nodecay_dom),
        ("test_domenico_lineardecay_model", testing_massbalance_lindecay_dom),
        ("test_domenico_instantreaction_model", testing_massbalance_instant_dom),
        ("test_domenico_instantreaction_model_inf", testing_massbalance_instant_dom_inf),
        (ModelParameters(), TypeError),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_balance_numerical(model, expected, request) -> None:
    """Test if mass balance is correctly calculated by comparing to precomputed results for Domenico model."""
    if isinstance(expected, dict):
        model_object = request.getfixturevalue(model)
        dictionary = mass_balance(model_object, time=3 * 365)
        for key, output_item in dictionary.items():
            assert expected[key] == pytest.approx(output_item)
    else:
        with pytest.raises(expected):
            mass_balance(model, time=3 * 365)
