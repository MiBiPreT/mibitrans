import pytest
from mibitrans.data.check_input import DomainValueError
from mibitrans.data.check_input import MissingValueError
from mibitrans.data.read import AttenuationParameters
from mibitrans.transport.domenico import InstantReaction
from mibitrans.transport.domenico import LinearDecay
from mibitrans.transport.domenico import NoDecay
from tests.test_example_data import test_att_pars
from tests.test_example_data import test_hydro_pars
from tests.test_example_data import test_model_pars
from tests.test_example_data import test_source_pars
from tests.test_example_data import testingdata_instantreaction
from tests.test_example_data import testingdata_lineardecay
from tests.test_example_data import testingdata_nodecay


@pytest.mark.parametrize(
    "att, error",
    [
        (AttenuationParameters(decay_rate=1), None),
        (AttenuationParameters(half_life=1), None),
        (
            AttenuationParameters(
                half_life=1, delta_oxygen=1.65, delta_nitrate=0.7, ferrous_iron=16.6, delta_sulfate=22.4, methane=6.6
            ),
            None,
        ),
    ],
)
def test_require_degradation_linear(att, error):
    """Test if LinearDecay class correctly raises error when correct degradation parameters are missing."""
    LinearDecay(test_hydro_pars, att, test_source_pars, test_model_pars)


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
    """Test if InstantReaction class correctly raises error when correct attenuation parameters are missing."""
    if error is None:
        InstantReaction(test_hydro_pars, att, test_source_pars, test_model_pars)
    else:
        with pytest.raises(error):
            InstantReaction(test_hydro_pars, att, test_source_pars, test_model_pars)


@pytest.mark.parametrize(
    "model, expected",
    [
        (NoDecay(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars), testingdata_nodecay),
        (
            LinearDecay(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars),
            testingdata_lineardecay,
        ),
        (
            InstantReaction(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars),
            testingdata_instantreaction,
        ),
    ],
)
def test_transport_equations_numerical(model, expected):
    """Test numerical output of transport equation child classes of Domenico."""
    assert model.cxyt == pytest.approx(expected)


@pytest.mark.parametrize(
    "x, y, t, expected",
    [
        (16, 0, 393, 2.83828613605873),
        (24, -5, 283, 0.5974811505254043),
        (-16, 0, 393, DomainValueError),
        ("nonsense", 0, 393, TypeError),
        (16, "nonsense", 393, TypeError),
        (16, 0, -10, DomainValueError),
        (16, 0, "nonsense", TypeError),
    ],
)
def test_domenico_sample(x, y, t, expected):
    """Tests if sample method from Domenico class works correctly."""
    model = NoDecay(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    if isinstance(expected, float):
        assert model.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            model.sample(x, y, t)
    elif isinstance(expected, list):
        with pytest.warns(expected[0]):
            assert model.sample(x, y, t) == pytest.approx(expected[1])
