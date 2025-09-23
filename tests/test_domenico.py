import numpy as np
import pytest
from mibitrans.data.read import AdsorptionParameters
from mibitrans.data.read import DegradationParameters
from mibitrans.data.read import ModelParameters
from mibitrans.transport.domenico import Domenico
from mibitrans.transport.domenico import InstantReaction
from mibitrans.transport.domenico import LinearDecay
from mibitrans.transport.domenico import NoDecay
from tests.test_example_data import test_ads_pars
from tests.test_example_data import test_deg_pars
from tests.test_example_data import test_hydro_pars
from tests.test_example_data import test_model_pars
from tests.test_example_data import test_source_pars
from tests.test_example_data import testingdata_instantreaction
from tests.test_example_data import testingdata_lineardecay
from tests.test_example_data import testingdata_nodecay

short_width_model_pars = ModelParameters(model_length=50, model_width=30, model_time=3 * 365, dx=10, dy=5, dt=1 * 365)
short_width_model_pars.model_width = test_source_pars.source_zone_boundary[-1] - 1


@pytest.mark.parametrize(
    "hydro, ads, source, model, error",
    [
        (test_hydro_pars, test_ads_pars, test_source_pars, test_model_pars, None),
        (1, test_ads_pars, test_source_pars, test_model_pars, TypeError),
        (test_hydro_pars, test_hydro_pars, test_source_pars, test_model_pars, TypeError),
        (test_hydro_pars, test_ads_pars, "wrong", test_model_pars, TypeError),
        (test_hydro_pars, test_ads_pars, test_source_pars, test_deg_pars, TypeError),
        (test_hydro_pars, test_ads_pars, test_source_pars, short_width_model_pars, UserWarning),
    ],
)
@pytest.mark.filterwarnings("ignore:UserWarning")
def test_domenico_parent(hydro, ads, source, model, error) -> None:
    """Test functionality, results and errors of Domenico parent class."""
    if error is None:
        parent = Domenico(hydro, ads, source, model)
        shape_arrays = (len(parent.t), len(parent.y), len(parent.x))
        # Source zone concentrations adapted for superposition should still have the same length as those in input
        assert (len(parent.c_source) == len(source.source_zone_concentration)) and (
            len(parent.c_source) == len(source.source_zone_boundary)
        )
        # Extent of y-domain should be at least the size of
        assert (np.max(parent.y) + abs(np.min(parent.y))) >= (np.max(source.source_zone_boundary) * 2)
        assert parent.xxx.shape == shape_arrays
        assert parent.yyy.shape == shape_arrays
        assert parent.ttt.shape == shape_arrays
        assert parent.ads_pars.retardation is not None
        assert hydro.velocity / parent.ads_pars.retardation == parent.rv
    elif error is UserWarning:
        with pytest.warns(UserWarning):
            parent = Domenico(hydro, ads, source, model)
            range_y = abs(parent.y[0]) + abs(parent.y[-1])
            assert range_y >= parent.src_pars.source_zone_boundary[-1] * 2
    elif error is TypeError:
        with pytest.raises(error):
            Domenico(hydro, ads, source, model)


@pytest.mark.parametrize(
    "ads, expected",
    [
        (AdsorptionParameters(retardation=1), 1),
        (AdsorptionParameters(bulk_density=2, partition_coefficient=10, fraction_organic_carbon=0.03), 3.4),
    ],
)
def test_retardation_calculation(ads, expected):
    """Test if retardation is calculated correctly when Domenico class is initialized."""
    parent = Domenico(test_hydro_pars, ads, test_source_pars, test_model_pars)
    assert parent.ads_pars.retardation == expected


@pytest.mark.parametrize(
    "deg, error",
    [
        (DegradationParameters(decay_rate=1), None),
        (DegradationParameters(half_life=1), None),
        (
            DegradationParameters(
                half_life=1, delta_oxygen=1.65, delta_nitrate=0.7, ferrous_iron=16.6, delta_sulfate=22.4, methane=6.6
            ),
            None,
        ),
        (
            DegradationParameters(
                delta_oxygen=1.65, delta_nitrate=0.7, ferrous_iron=16.6, delta_sulfate=22.4, methane=6.6
            ),
            ValueError,
        ),
    ],
)
def test_require_degradation_linear(deg, error):
    """Test if LinearDecay class correctly raises error when correct degradation parameters are missing."""
    if error is None:
        LinearDecay(test_hydro_pars, test_ads_pars, deg, test_source_pars, test_model_pars)
    else:
        with pytest.raises(error):
            LinearDecay(test_hydro_pars, test_ads_pars, deg, test_source_pars, test_model_pars)


@pytest.mark.parametrize(
    "deg, error",
    [
        (
            DegradationParameters(
                half_life=1, delta_oxygen=1.65, delta_nitrate=0.7, ferrous_iron=16.6, delta_sulfate=22.4, methane=6.6
            ),
            None,
        ),
        (
            DegradationParameters(
                delta_oxygen=1.65, delta_nitrate=0.7, ferrous_iron=16.6, delta_sulfate=22.4, methane=6.6
            ),
            None,
        ),
        (DegradationParameters(decay_rate=1), ValueError),
        (DegradationParameters(half_life=1), ValueError),
        (
            DegradationParameters(half_life=1, delta_oxygen=1.65, delta_nitrate=0.7, ferrous_iron=16.6, methane=6.6),
            ValueError,
        ),
    ],
)
def test_require_degradation_instant(deg, error):
    """Test if InstantReaction class correctly raises error when correct degradation parameters are missing."""
    if error is None:
        InstantReaction(test_hydro_pars, test_ads_pars, deg, test_source_pars, test_model_pars)
    else:
        with pytest.raises(error):
            InstantReaction(test_hydro_pars, test_ads_pars, deg, test_source_pars, test_model_pars)


@pytest.mark.parametrize(
    "model, expected",
    [
        (NoDecay(test_hydro_pars, test_ads_pars, test_source_pars, test_model_pars), testingdata_nodecay),
        (
            LinearDecay(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars),
            testingdata_lineardecay,
        ),
        (
            InstantReaction(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars),
            testingdata_instantreaction,
        ),
    ],
)
def test_transport_equations_numerical(model, expected):
    """Test numerical output of transport equation child classes of Domenico."""
    assert model.cxyt == pytest.approx(expected)
