import pytest
from mibitrans.data.check_input import DomainValueError
from mibitrans.data.check_input import MissingValueError
from mibitrans.data.read import AttenuationParameters
from mibitrans.transport.karanovic import Instant
from mibitrans.transport.karanovic import Linear
from tests.test_example_data import test_att_pars
from tests.test_example_data import test_hydro_pars
from tests.test_example_data import test_model_pars
from tests.test_example_data import test_source_pars
from tests.test_example_data import testingdata_instantreaction_karanovic
from tests.test_example_data import testingdata_lineardecay_karanovic


@pytest.mark.parametrize(
    "att, error",
    [
        (
            AttenuationParameters(
                half_life=1, delta_oxygen=1.65, delta_nitrate=0.7, ferrous_iron=16.6, delta_sulfate=22.4, methane=6.6
            ),
            None,
        ),
        (
            AttenuationParameters(
                delta_oxygen=1.65, delta_nitrate=0.7, ferrous_iron=16.6, delta_sulfate=22.4, methane=6.6
            ),
            None,
        ),
        (AttenuationParameters(decay_rate=1), MissingValueError),
        (AttenuationParameters(half_life=1), MissingValueError),
        (
            AttenuationParameters(half_life=1, delta_oxygen=1.65, delta_nitrate=0.7, ferrous_iron=16.6, methane=6.6),
            MissingValueError,
        ),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_require_degradation_instant(att, error):
    """Test if Instan class correctly raises error when correct attenuation parameters are missing."""
    if error is None:
        Instant(test_hydro_pars, att, test_source_pars, test_model_pars)
    else:
        with pytest.raises(error):
            Instant(test_hydro_pars, att, test_source_pars, test_model_pars)


@pytest.mark.parametrize(
    "model, expected",
    [
        (
            Linear(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars),
            testingdata_lineardecay_karanovic,
        ),
        (
            Instant(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars),
            testingdata_instantreaction_karanovic,
        ),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_transport_equations_numerical(model, expected):
    """Test numerical output of transport equation child classes of Karanovic."""
    assert model.cxyt == pytest.approx(expected)


sample_model = Linear(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)


@pytest.mark.parametrize(
    "x, y, t, expected",
    [
        (16, 0, 393, 0.2980583920684923),
        (24, -5, 283, 0.03725197246248769),
        (-16, 0, 393, DomainValueError),
        ("nonsense", 0, 393, TypeError),
        (16, "nonsense", 393, TypeError),
        (16, 0, -10, DomainValueError),
        (16, 0, "nonsense", TypeError),
    ],
)
def test_karanovic_sample(x, y, t, expected):
    """Tests if sample method from Karanovic class works correctly."""
    if isinstance(expected, float):
        assert sample_model.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            sample_model.sample(x, y, t)
    elif isinstance(expected, list):
        with pytest.warns(expected[0]):
            assert sample_model.sample(x, y, t) == pytest.approx(expected[1])
