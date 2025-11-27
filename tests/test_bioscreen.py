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
def test_transport_equations_numerical(model, expected, request):
    """Test numerical output of transport equation child classes of Domenico."""
    model_object = request.getfixturevalue(model)
    assert model_object.cxyt == pytest.approx(expected)


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
def test_anatrans_sample_linear(x, y, t, expected, test_bioscreen_model_lineardecay):
    """Tests if sample method from Domenico class works correctly."""
    if isinstance(expected, float):
        assert test_bioscreen_model_lineardecay.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            test_bioscreen_model_lineardecay.sample(x, y, t)


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
def test_anatrans_sample_instantreaction(x, y, t, expected, test_bioscreen_model_instantreaction):
    """Tests if sample method from Domenico class works correctly."""
    if isinstance(expected, float):
        assert test_bioscreen_model_instantreaction.sample(x, y, t) == pytest.approx(expected)
    elif expected is ValueError or expected is TypeError:
        with pytest.raises(expected):
            test_bioscreen_model_instantreaction.sample(x, y, t)
