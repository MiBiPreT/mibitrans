import pytest
from mibitrans.data.check_input import DomainValueError
from tests.test_example_data import testingdata_instantreaction_bioscreen
from tests.test_example_data import testingdata_lineardecay_bioscreen
from tests.test_example_data import testingdata_nodecay_bioscreen


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_bioscreen_model_nodecay", testingdata_nodecay_bioscreen),
        ("test_bioscreen_model_lineardecay", testingdata_lineardecay_bioscreen),
        ("test_bioscreen_model_instantreaction", testingdata_instantreaction_bioscreen),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_transport_equation_numerical_bioscreen(model, expected, request):
    """Test numerical output of transport equation of Anatrans, by comparing to pre-calculated values."""
    model, results = request.getfixturevalue(model)
    assert model.cxyt == pytest.approx(expected)
    assert results.cxyt == pytest.approx(expected)


@pytest.mark.parametrize(
    "x, y, t, expected",
    [
        (9, 0, 629, 1.2260728205395477),
        (15, -7, 256, 0.21033402922523056),
        (-16, 0, 393, DomainValueError),
        ("nonsense", 0, 393, TypeError),
        (16, "nonsense", 393, TypeError),
        (16, 0, -10, DomainValueError),
        (16, 0, "nonsense", TypeError),
    ],
)
def test_bioscreen_sample_linear(x, y, t, expected, test_bioscreen_model_lineardecay):
    """Test if sample method from Bioscreen works correctly, and gives expected output for linear models."""
    model, results = test_bioscreen_model_lineardecay
    if isinstance(expected, float):
        assert model.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            model.sample(x, y, t)


@pytest.mark.parametrize(
    "x, y, t, expected",
    [
        (13, 0, 354, 2.721953070355462),
        (11, 3, 752, 5.01432266465888),
        (-16, 0, 393, DomainValueError),
        ("nonsense", 0, 393, TypeError),
        (16, "nonsense", 393, TypeError),
        (16, 0, -10, DomainValueError),
        (16, 0, "nonsense", TypeError),
    ],
)
def test_bioscreen_sample_instantreaction(x, y, t, expected, test_bioscreen_model_instantreaction):
    """Test if sample method from Bioscreen works correctly, and gives expected output for instant reaction models."""
    model, results = test_bioscreen_model_instantreaction
    if isinstance(expected, float):
        assert model.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            model.sample(x, y, t)
