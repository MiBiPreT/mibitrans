import pytest
from mibitrans.data.check_input import DomainValueError
from mibitrans.data.parameter_information import InstantElectronAcceptors
from mibitrans.data.parameter_information import UtilizationFactor
from mibitrans.transport.model_parent import Results
from mibitrans.transport.models import Bioscreen


@pytest.mark.parametrize(
    "pars, expected",
    [
        ({"electron_acceptors": [0.2, 0.4, 1, 0.5, 1]}, None),
        ({"electron_acceptors": [0.2, 0.4, 1, 0.5, 1], "utilization_factor": [2.1, 1, 2, 3, 0.2]}, None),
        # Test for accepting ElectronAcceptors and UtilizationFactor dataclasses
        (
            {
                "electron_acceptors": InstantElectronAcceptors(0.2, 0.4, 1, 0.5, 1),
                "utilization_factor": UtilizationFactor(2.1, 1, 2, 3, 0.2),
            },
            None,
        ),
        # Test for accepting dictionaries
        (
            {
                "electron_acceptors": dict(
                    delta_oxygen=0.2, delta_nitrate=0.4, ferrous_iron=1, delta_sulfate=0.5, methane=1
                ),
                "utilization_factor": dict(
                    util_oxygen=2.1, util_nitrate=1, util_ferrous_iron=2, util_sulfate=3, util_methane=0.2
                ),
            },
            None,
        ),
        ({"utilization_factor": [2.1, 1, 2, 3, 0.2]}, TypeError),
        ({"electron_acceptors": "No ea for you", "utilization_factor": [2.1, 1, 2, 3, 0.2]}, TypeError),
        ({"electron_acceptors": [0.2, 0.4, 1, 0.5, 1], "utilization_factor": 3}, TypeError),
        ({"electron_acceptors": [0.2, "horseradish", 1, 0.5, 1], "utilization_factor": [2.1, 1, 2, 3, 0.2]}, TypeError),
        ({"electron_acceptors": [0.2, 0.4, 1, 0.5, 1], "utilization_factor": [2.1, 1, -2, 3, 0.2]}, DomainValueError),
    ],
)
def test_instant_reaction_setup(
    pars, expected, test_hydro_pars, test_att_pars, test_source_pars, test_model_pars
) -> None:
    """Test if the instant reaction method accepts correct data and raises the intended errors if not."""
    model_object = Bioscreen(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    if expected is None:
        model_object.instant_reaction(**pars)
        assert model_object.mode == "instant_reaction", (
            f"Model mode should be 'instant_reaction', but is '{model_object.mode}' instead."
        )
        assert model_object.biodegradation_capacity is not None, "No biodegradation capacity was calculated."
    elif expected in [TypeError, ValueError]:
        with pytest.raises(expected):
            model_object.instant_reaction(**pars)


@pytest.mark.parametrize(
    "model, expected",
    [
        ("test_bioscreen_model_nodecay", "nodecay_bioscreen"),
        ("test_bioscreen_model_lineardecay", "lineardecay_bioscreen"),
        ("test_bioscreen_model_instantreaction", "instantreaction_bioscreen"),
        ("test_bioscreen_model_fringe", "instantreaction_bioscreen"),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_transport_equation_numerical_bioscreen(model, expected, request, test_example_data):
    """Test numerical output of transport equation of Anatrans, by comparing to pre-calculated values."""
    mod, results = request.getfixturevalue(model)
    expected_value = getattr(test_example_data, expected)
    if model == "test_bioscreen_model_fringe":
        assert mod.cxyt[0] == pytest.approx(expected_value)
        assert results.cxyt[0] == pytest.approx(expected_value)
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
