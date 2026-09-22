import pytest
from mibitrans.data.check_input import DomainValueError
from mibitrans.transport.model_parent import Results
from mibitrans.transport.models import Mibitrans


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_mibitrans_model_nodecay", "nodecay_mibitrans"),
        ("test_mibitrans_model_lineardecay", "lineardecay_mibitrans"),
        ("test_mibitrans_model_instantreaction", "instantreaction_mibitrans"),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_transport_equation_numerical_mibitrans(model, expected, request, test_example_data):
    """Test numerical output of transport equation of Mibitrans, by comparing to pre-calculated values."""
    model, results = request.getfixturevalue(model)
    expected_value = getattr(test_example_data, expected)
    assert model.cxyt == pytest.approx(expected_value)
    assert results.cxyt == pytest.approx(expected_value)


def test_transport_chain_decay_runs(test_hydro_pars, test_att_pars_chain, test_source_pars_chain, test_model_pars):
    """Test if running the chain decay method produces expected Results object."""
    model_obj = Mibitrans(test_hydro_pars, test_att_pars_chain, test_source_pars_chain, test_model_pars)
    model_obj.chain_decay([0.8, 0.7])
    results = model_obj.run()
    assert isinstance(results, Results), "Result object is not of type Results"


@pytest.mark.parametrize(
    "x, y, t, expected",
    [
        (16, 0, 393, 4.6652115692165195),
        (24, -5, 527, 1.7909410090260753),
        (-16, 0, 393, DomainValueError),
        ("nonsense", 0, 393, TypeError),
        (16, "nonsense", 393, TypeError),
        (16, 0, -10, DomainValueError),
        (16, 0, "nonsense", TypeError),
    ],
)
def test_mibitrans_linear_sample(x, y, t, expected, test_mibitrans_model_lineardecay):
    """Test if sample method from Mibitrans works correctly, and gives expected output for linear models."""
    model, results = test_mibitrans_model_lineardecay
    if isinstance(expected, float):
        assert model.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            model.sample(x, y, t)


@pytest.mark.parametrize(
    "x, y, t, expected",
    [
        (20, 0, 476, 6.3213083960634435),
        (35, 7, 745, 2.5356628358944633),
        (-16, 0, 393, DomainValueError),
        ("nonsense", 0, 393, TypeError),
        (16, "nonsense", 393, TypeError),
        (16, 0, -10, DomainValueError),
        (16, 0, "nonsense", TypeError),
    ],
)
def test_mibitrans_instant_sample(x, y, t, expected, test_mibitrans_model_instantreaction):
    """Test if sample method from Mibitrans works correctly, and gives expected output for instant reaction models."""
    model, results = test_mibitrans_model_instantreaction
    if isinstance(expected, float):
        assert model.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            model.sample(x, y, t)
