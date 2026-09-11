import pytest
from mibitrans.data.check_input import DomainValueError
from mibitrans.transport.model_parent import Results
from mibitrans.transport.models import Anatrans


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_anatrans_model_nodecay", "nodecay_anatrans"),
        ("test_anatrans_model_lineardecay", "lineardecay_anatrans"),
        ("test_anatrans_model_instantreaction", "instantreaction_anatrans"),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_transport_equation_numerical_anatrans(model, expected, request, test_example_data):
    """Test numerical output of transport equation of Anatrans, by comparing to pre-calculated values."""
    model, results = request.getfixturevalue(model)
    expected_value = getattr(test_example_data, expected)
    assert model.cxyt == pytest.approx(expected_value)
    assert results.cxyt == pytest.approx(expected_value)


def test_transport_chain_decay_runs(test_hydro_pars, test_att_pars_chain, test_source_pars_chain, test_model_pars):
    """Test if running the chain decay method produces expected Results object."""
    model_obj = Anatrans(test_hydro_pars, test_att_pars_chain, test_source_pars_chain, test_model_pars)
    model_obj.chain_decay([0.8, 0.7])
    results = model_obj.run()
    assert isinstance(results, Results), "Result object is not of type Results"


@pytest.mark.parametrize(
    "x, y, t, expected",
    [
        (16, 0, 393, 4.041372051306399),
        (24, -5, 283, 1.5760137786262713),
        (-16, 0, 393, DomainValueError),
        ("nonsense", 0, 393, TypeError),
        (16, "nonsense", 393, TypeError),
        (16, 0, -10, DomainValueError),
        (16, 0, "nonsense", TypeError),
    ],
)
def test_anatrans_sample_linear(x, y, t, expected, test_anatrans_model_lineardecay):
    """Test if sample method from Anatrans works correctly, and gives expected output for linear models."""
    model, results = test_anatrans_model_lineardecay
    if isinstance(expected, float):
        assert model.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            model.sample(x, y, t)


@pytest.mark.parametrize(
    "x, y, t, expected",
    [
        (20, 0, 476, 5.540354132380653),
        (54, 3, 1045, 3.501165953555688),
        (-16, 0, 393, DomainValueError),
        ("nonsense", 0, 393, TypeError),
        (16, "nonsense", 393, TypeError),
        (16, 0, -10, DomainValueError),
        (16, 0, "nonsense", TypeError),
    ],
)
def test_anatrans_sample_instantreaction(x, y, t, expected, test_anatrans_model_instantreaction):
    """Test if sample method from Anatrans works correctly, and gives expected output for instant reaction models."""
    model, results = test_anatrans_model_instantreaction
    if isinstance(expected, float):
        assert model.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            model.sample(x, y, t)
