"""Author: Jorrit Bakker.

File testing functionality of mass_balance module.
"""

import numpy as np
import pytest
from mibitrans.data.parameters import ModelParameters
from mibitrans.data.parameters import SourceParameters
from mibitrans.transport.models import Anatrans
from mibitrans.transport.models import Bioscreen
from mibitrans.transport.models import Mibitrans
from tests.test_example_data import testing_massbalance_instant_ana_inf
from tests.test_example_data import testing_massbalance_instant_mbt
from tests.test_example_data import testing_massbalance_instant_mbt_inf
from tests.test_example_data import testing_massbalance_lineardecay_ana
from tests.test_example_data import testing_massbalance_nodecay_bio


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
def test_bioscreen_nodecay_model_mb(test_hydro_pars, test_att_pars_nodecay, test_source_pars, test_model_pars):
    """Bioscreen with linear decay fixture mass balance object for testing."""
    obj = Bioscreen(test_hydro_pars, test_att_pars_nodecay, test_source_pars, test_model_pars)
    res = obj.run()
    return res.mass_balance()


@pytest.fixture(scope="module")
def test_anatrans_lineardecay_model_mb(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """Anatrans with linear decay fixture mass balance object for testing."""
    obj = Anatrans(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    res = obj.run()
    return res.mass_balance()


@pytest.fixture(scope="module")
def test_mibitrans_instantreaction_model_mb(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """Mibitrans with instant reaction fixture mass balance object for testing."""
    obj = Mibitrans(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    obj.instant_reaction(dict(delta_oxygen=0.5, delta_nitrate=0.5, ferrous_iron=0.5, delta_sulfate=0.5, methane=0.5))
    res = obj.run()
    return res.mass_balance()


@pytest.fixture(scope="module")
def test_mibitrans_instantreaction_model_mb_inf(test_hydro_pars, test_att_pars, test_source_pars_inf, test_model_pars):
    """Mibitrans with instant reaction and infinite source mass fixture mass balance  object for testing."""
    obj = Mibitrans(test_hydro_pars, test_att_pars, test_source_pars_inf, test_model_pars)
    obj.instant_reaction(dict(delta_oxygen=0.5, delta_nitrate=0.5, ferrous_iron=0.5, delta_sulfate=0.5, methane=0.5))
    res = obj.run()
    return res.mass_balance()


@pytest.fixture(scope="module")
def test_anatrans_instantreaction_model_mb_inf(test_hydro_pars, test_att_pars, test_source_pars_inf, test_model_pars):
    """Anatrans with instant reaction and infinite source mass fixture mass balance  object for testing."""
    obj = Anatrans(test_hydro_pars, test_att_pars, test_source_pars_inf, test_model_pars)
    obj.instant_reaction(dict(delta_oxygen=0.5, delta_nitrate=0.5, ferrous_iron=0.5, delta_sulfate=0.5, methane=0.5))
    res = obj.run()
    return res.mass_balance()


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_bioscreen_nodecay_model_mb", testing_massbalance_nodecay_bio),
        ("test_anatrans_lineardecay_model_mb", testing_massbalance_lineardecay_ana),
        ("test_mibitrans_instantreaction_model_mb", testing_massbalance_instant_mbt),
        ("test_mibitrans_instantreaction_model_mb_inf", testing_massbalance_instant_mbt_inf),
        ("test_anatrans_instantreaction_model_mb_inf", testing_massbalance_instant_ana_inf),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
@pytest.mark.filterwarnings("ignore:Contaminant plume extents")
def test_balance_numerical_mibitrans(model, expected, request) -> None:
    """Test if mass balance is correctly calculated by comparing to precomputed results for Mibitrans model."""
    mb_object = request.getfixturevalue(model)
    dictionary = mb_object.__dict__
    for key, output_item in expected.items():
        if isinstance(output_item, dict):
            for key_2, value_2 in output_item.items():
                assert value_2 == pytest.approx(dictionary[key][key_2])
        else:
            assert output_item == pytest.approx(dictionary[key]), f"mb_item {key}"

    mb_object(2 * 365)
    dictionary_single = mb_object.__dict__
    for key, output_item in expected.items():
        if isinstance(output_item, np.ndarray):
            assert output_item[1] == pytest.approx(dictionary_single[key]), f"mb_item {key}"
        elif isinstance(output_item, dict):
            for key_2, value_2 in output_item.items():
                assert value_2[1] == pytest.approx(dictionary_single[key][key_2]), f"mb_item {key}"
        else:
            assert output_item == pytest.approx(dictionary_single[key]), f"mb_item {key}"
