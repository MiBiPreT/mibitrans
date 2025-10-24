import pytest
from mibitrans.data.check_input import DomainValueError
from mibitrans.data.check_input import MissingValueError
from mibitrans.data.parameters import AttenuationParameters
from mibitrans.transport.karanovic import InstantReaction
from tests.test_example_data import testingdata_instantreaction_karanovic
from tests.test_example_data import testingdata_lineardecay_karanovic
from tests.test_example_data import testingdata_nodecay_karanovic


@pytest.mark.parametrize(
    "att, error",
    [
        (
            AttenuationParameters(
                half_life=1,
                electron_acceptors=dict(
                    delta_oxygen=1.65, delta_nitrate=0.7, ferrous_iron=16.6, delta_sulfate=22.4, methane=6.6
                ),
            ),
            None,
        ),
        (
            AttenuationParameters(
                electron_acceptors=dict(
                    delta_oxygen=1.65, delta_nitrate=0.7, ferrous_iron=16.6, delta_sulfate=22.4, methane=6.6
                )
            ),
            None,
        ),
        (AttenuationParameters(decay_rate=1), MissingValueError),
        (AttenuationParameters(half_life=1), MissingValueError),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_require_degradation_instant(att, error, test_hydro_pars, test_source_pars, test_model_pars):
    """Test if Instan class correctly raises error when correct attenuation parameters are missing."""
    if error is None:
        InstantReaction(test_hydro_pars, att, test_source_pars, test_model_pars)
    else:
        with pytest.raises(error):
            InstantReaction(test_hydro_pars, att, test_source_pars, test_model_pars)


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_karanovic_nodecay_model", testingdata_nodecay_karanovic),
        ("test_karanovic_lineardecay_model", testingdata_lineardecay_karanovic),
        ("test_karanovic_instantreaction_model", testingdata_instantreaction_karanovic),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_transport_equations_numerical(model, expected, request):
    """Test numerical output of transport equation child classes of Karanovic."""
    model_object = request.getfixturevalue(model)
    assert model_object.cxyt == pytest.approx(expected)


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
def test_karanovic_linear_sample(x, y, t, expected, test_karanovic_lineardecay_model):
    """Tests if sample method from Karanovic class works correctly for LinearDecay."""
    if isinstance(expected, float):
        assert test_karanovic_lineardecay_model.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            test_karanovic_lineardecay_model.sample(x, y, t)
    elif isinstance(expected, list):
        with pytest.warns(expected[0]):
            assert test_karanovic_lineardecay_model.sample(x, y, t) == pytest.approx(expected[1])


@pytest.mark.parametrize(
    "x, y, t, expected",
    [
        (20, 0, 476, 3.8101869779573443),
        (11, 7, 193, 2.0276832492832924),
        (-16, 0, 393, DomainValueError),
        ("nonsense", 0, 393, TypeError),
        (16, "nonsense", 393, TypeError),
        (16, 0, -10, DomainValueError),
        (16, 0, "nonsense", TypeError),
    ],
)
def test_karanovic_instant_sample(x, y, t, expected, test_karanovic_instantreaction_model):
    """Tests if sample method from Karanovic class works correctly for InstantReaction."""
    if isinstance(expected, float):
        assert test_karanovic_instantreaction_model.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            test_karanovic_instantreaction_model.sample(x, y, t)
    elif isinstance(expected, list):
        with pytest.warns(expected[0]):
            assert test_karanovic_instantreaction_model.sample(x, y, t) == pytest.approx(expected[1])
