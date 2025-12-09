import pytest
from mibitrans.data.check_input import DomainValueError
from tests.test_example_data import testingdata_instantreaction_anatrans
from tests.test_example_data import testingdata_lineardecay_anatrans
from tests.test_example_data import testingdata_nodecay_anatrans


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_anatrans_model_nodecay", testingdata_nodecay_anatrans),
        ("test_anatrans_model_lineardecay", testingdata_lineardecay_anatrans),
        ("test_anatrans_model_instantreaction", testingdata_instantreaction_anatrans),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_transport_equation_numerical_anatrans(model, expected, request):
    """Test numerical output of transport equation of Anatrans, by comparing to pre-calculated values."""
    model, results = request.getfixturevalue(model)
    assert model.cxyt == pytest.approx(expected)
    assert results.cxyt == pytest.approx(expected)


@pytest.mark.parametrize(
    "x, y, t, expected",
    [
        (16, 0, 393, 0.2403406438598838),
        (24, -5, 283, 0.031529981399875194),
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
        (20, 0, 476, 3.076798202181921),
        (11, 7, 193, 2.0747279062256183),
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
