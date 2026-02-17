import pytest
from mibitrans.data.check_input import DomainValueError
from tests.test_example_data import testingdata_instantreaction_mibitrans
from tests.test_example_data import testingdata_lineardecay_mibitrans
from tests.test_example_data import testingdata_nodecay_mibitrans


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_mibitrans_model_nodecay", testingdata_nodecay_mibitrans),
        ("test_mibitrans_model_lineardecay", testingdata_lineardecay_mibitrans),
        ("test_mibitrans_model_instantreaction", testingdata_instantreaction_mibitrans),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_transport_equation_numerical_mibitrans(model, expected, request):
    """Test numerical output of transport equation of Mibitrans, by comparing to pre-calculated values."""
    model, results = request.getfixturevalue(model)
    assert model.cxyt == pytest.approx(expected)
    assert results.cxyt == pytest.approx(expected)


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
