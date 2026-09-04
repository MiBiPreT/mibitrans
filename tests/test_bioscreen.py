import pytest
from mibitrans.data.check_input import DomainValueError
from mibitrans.transport.model_parent import Results
from mibitrans.transport.models import Bioscreen


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_bioscreen_model_nodecay", "nodecay_bioscreen"),
        ("test_bioscreen_model_lineardecay", "lineardecay_bioscreen"),
        ("test_bioscreen_model_instantreaction", "instantreaction_bioscreen"),
        ("test_bioscreen_model_fringe", "instantreaction_bioscreen")
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_transport_equation_numerical_bioscreen(model, expected, request, test_example_data):
    """Test numerical output of transport equation of Anatrans, by comparing to pre-calculated values."""
    mod, results = request.getfixturevalue(model)
    expected_value = getattr(test_example_data, expected)
    if model == "test_bioscreen_model_fringe":
        # Until decision on source depletion fringe degradation solution, skip testing output
        pass
        # assert mod.cxyt[0] == pytest.approx(expected_value)
        # assert results.cxyt[0] == pytest.approx(expected_value)
    else:
        assert mod.cxyt == pytest.approx(expected_value)
        assert results.cxyt == pytest.approx(expected_value)


def test_transport_chain_decay_runs(test_hydro_pars, test_att_pars_chain, test_source_pars_chain, test_model_pars):
    """Test if running the chain decay method produces expected Results object."""
    model_obj = Bioscreen(test_hydro_pars, test_att_pars_chain, test_source_pars_chain, test_model_pars)
    model_obj.chain_decay([0.8, 0.7])
    results = model_obj.run()
    assert isinstance(results, Results), "Result object is not of type Results"


@pytest.mark.parametrize(
    "x, y, t, expected",
    [
        (9, 0, 629, 6.222919410416837),
        (15, -7, 256, 1.5292214426149926),
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
        (13, 0, 354, 5.134230309076454),
        (11, 3, 752, 5.645872356690198),
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
