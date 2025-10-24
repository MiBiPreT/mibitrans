import numpy as np
import pytest
from mibitrans.data.check_input import DomainValueError
from mibitrans.data.check_input import MissingValueError
from mibitrans.data.parameters import AttenuationParameters
from mibitrans.transport import domenico as dom
from tests.test_example_data import testingdata_instantreaction
from tests.test_example_data import testingdata_lineardecay
from tests.test_example_data import testingdata_nodecay


@pytest.mark.parametrize(
    "model_type",
    [
        "no_decay",
        "linear",
        "instant",
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_auto_run_off_and_reset(model_type, test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """Test if model does not auto run when set to False and resets whenever property gets changed."""
    match model_type:
        case "no_decay":
            mod = dom.NoDecay(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars, auto_run=False)
        case "linear":
            mod = dom.LinearDecay(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars, auto_run=False)
        case "instant":
            mod = dom.InstantReaction(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars, auto_run=False)

    assert not mod.has_run
    assert np.sum(mod.cxyt) == 0.0

    mod.run()

    assert mod.has_run
    assert np.sum(mod.cxyt) > 0.0

    mod.hydrological_parameters.velocity = 10

    assert not mod.has_run
    assert np.sum(mod.cxyt) == 0.0


@pytest.mark.parametrize(
    "att, error",
    [
        (AttenuationParameters(decay_rate=1), None),
        (AttenuationParameters(half_life=1), None),
        (
            AttenuationParameters(
                half_life=1,
                electron_acceptors=dict(
                    delta_oxygen=1.65, delta_nitrate=0.7, ferrous_iron=16.6, delta_sulfate=22.4, methane=6.6
                ),
            ),
            None,
        ),
    ],
)
def test_require_degradation_linear(att, error, test_hydro_pars, test_source_pars, test_model_pars):
    """Test if LinearDecay class correctly raises error when correct degradation parameters are missing."""
    dom.LinearDecay(test_hydro_pars, att, test_source_pars, test_model_pars)


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
    """Test if InstantReaction class correctly raises error when correct attenuation parameters are missing."""
    if error is None:
        dom.InstantReaction(test_hydro_pars, att, test_source_pars, test_model_pars)
    else:
        with pytest.raises(error):
            dom.InstantReaction(test_hydro_pars, att, test_source_pars, test_model_pars)


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_domenico_nodecay_model", testingdata_nodecay),
        ("test_domenico_lineardecay_model", testingdata_lineardecay),
        ("test_domenico_instantreaction_model", testingdata_instantreaction),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_transport_equations_numerical(model, expected, request):
    """Test numerical output of transport equation child classes of Domenico."""
    model_object = request.getfixturevalue(model)
    assert model_object.cxyt == pytest.approx(expected)


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
def test_domenico_sample(x, y, t, expected, test_domenico_nodecay_model):
    """Tests if sample method from Domenico class works correctly."""
    if isinstance(expected, float):
        assert test_domenico_nodecay_model.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            test_domenico_nodecay_model.sample(x, y, t)
    elif isinstance(expected, list):
        with pytest.warns(expected[0]):
            assert test_domenico_nodecay_model.sample(x, y, t) == pytest.approx(expected[1])
