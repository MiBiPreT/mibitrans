import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from mibitrans.data.check_input import DomainValueError
from mibitrans.data.parameter_information import ElectronAcceptors
from mibitrans.data.parameter_information import UtilizationFactor
from mibitrans.data.parameters import AttenuationParameters
from mibitrans.transport.model_parent import Transport3D


class Transport3DConcrete(Transport3D):
    """Concrete class of abstract class Transport3D, for purpose of testing."""

    def short_description(self):
        """Short string describing model type."""
        return "Test version of abstract class Transport3D"

    def run(self):
        """Method that runs the model and ensures that initialisation is performed."""
        self._calculate_concentration_for_all_xyt()

    def sample(self, x_position, y_position, t_position):
        """Method that calculates concentration at single, specified location in model domain."""
        return self.cxyt

    def _calculate_concentration_for_all_xyt(self) -> np.ndarray:
        """Method that calculates and return concentration array for all model x, y and t."""
        return self.cxyt


@pytest.mark.parametrize(
    "hydro, att, source, model, error",
    [
        ("test_hydro_pars", "test_att_pars", "test_source_pars", "test_model_pars", None),
        (1, "test_att_pars", "test_source_pars", "test_model_pars", TypeError),
        ("test_hydro_pars", "test_hydro_pars", "test_source_pars", "test_model_pars", TypeError),
        ("test_hydro_pars", "test_att_pars", dict(), "test_model_pars", TypeError),
        ("test_hydro_pars", "test_att_pars", "test_source_pars", "test_model_pars_short", UserWarning),
    ],
)
def test_transport_parent(hydro, att, source, model, error, request) -> None:
    """Test functionality, results and errors of Transport3D parent class."""
    args = []
    for entry in [hydro, att, source, model]:
        if isinstance(entry, str):
            args.append(request.getfixturevalue(entry))
        else:
            args.append(entry)

    if error is None:
        parent = Transport3DConcrete(*args)
        # Source zone concentrations adapted for superposition should still have the same length as those in input
        assert (len(parent.c_source) == len(args[2].source_zone_concentration)) and (
            len(parent.c_source) == len(args[2].source_zone_boundary)
        )
        # Extent of y-domain should be at least the size of
        assert (np.max(parent.y) + abs(np.min(parent.y))) >= (np.max(args[2].source_zone_boundary) * 2)
        assert parent.xxx.shape == (1, 1, len(parent.x))
        assert parent.yyy.shape == (1, len(parent.y), 1)
        assert parent.ttt.shape == (len(parent.t), 1, 1)
        assert parent._att_pars.retardation is not None
        assert args[0].velocity / parent._att_pars.retardation == parent.rv
    elif error is UserWarning:
        with pytest.warns(UserWarning):
            parent = Transport3DConcrete(*args)
            range_y = abs(parent.y[0]) + abs(parent.y[-1])
            assert range_y >= parent._src_pars.source_zone_boundary[-1] * 2
    elif error is TypeError:
        with pytest.raises(error):
            Transport3DConcrete(*args)


@pytest.mark.parametrize(
    "att, expected",
    [
        (AttenuationParameters(retardation=1), 1),
        (AttenuationParameters(bulk_density=2, partition_coefficient=10, fraction_organic_carbon=0.03), 3.4),
    ],
)
def test_retardation_calculation(att, expected, test_hydro_pars, test_source_pars, test_model_pars) -> None:
    """Test if retardation is calculated correctly when Transport3D class is initialized."""
    parent = Transport3DConcrete(test_hydro_pars, att, test_source_pars, test_model_pars)
    assert parent._att_pars.retardation == expected


@pytest.mark.parametrize(
    "pars, expected",
    [
        ({"electron_acceptors": [0.2, 0.4, 1, 0.5, 1]}, None),
        ({"electron_acceptors": [0.2, 0.4, 1, 0.5, 1], "utilization_factor": [2.1, 1, 2, 3, 0.2]}, None),
        # Test for accepting ElectronAcceptors and UtilizationFactor dataclasses
        (
            {
                "electron_acceptors": ElectronAcceptors(0.2, 0.4, 1, 0.5, 1),
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
    model_object = Transport3DConcrete(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    if expected is None:
        model_object.instant_reaction(**pars)
        assert model_object.mode == "instant_reaction", (
            f"Model mode should be 'instant_reaction', but is '{model_object.mode}' instead."
        )
        assert model_object.biodegradation_capacity is not None, "No biodegradation capacity was calculated."
    elif expected in [TypeError, ValueError]:
        with pytest.raises(expected):
            model_object.instant_reaction(**pars)


def test_mode_switch(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars) -> None:
    """Test if model correctly switches modes when using the mode property."""
    pars = {"electron_acceptors": [0.2, 0.4, 1, 0.5, 1], "utilization_factor": [2.1, 1, 2, 3, 0.2]}
    model_object = Transport3DConcrete(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    with pytest.raises(ValueError):
        # Model should raise error if attempting to switch to instant reaction model without providing parameters.
        model_object.mode = "instant_reaction"
    model_object.instant_reaction(**pars)
    assert (
        model_object.c_source[-1]
        == test_source_pars.source_zone_concentration[-1] + model_object.biodegradation_capacity
    ), ("Biodegradation capacity was not added to the outermost source zone.",)
    model_object.mode = "linear"
    assert model_object.mode == "linear", "Model mode was not switched to 'linear'."

def test_model_results_independance(test_mibitrans_model_instantreaction):
    """Test to make sure that changing model parameters does not change result parameters."""
    model, results = test_mibitrans_model_instantreaction
    assert model.hydrological_parameters.velocity == results.hydrological_parameters.velocity
    model.hydrological_parameters.velocity += 1
    assert model.hydrological_parameters.velocity != results.hydrological_parameters.velocity

def test_plotting_methods(test_anatrans_model_nodecay):
    """Test if plotting methods defined in parent model are working."""
    model, results = test_anatrans_model_nodecay

    results.centerline()
    assert isinstance(plt.gca(), matplotlib.axes._axes.Axes), "No or incorrect plot was created for centerline."
    plt.clf()

    results.transverse(x_position=1)
    assert isinstance(plt.gca(), matplotlib.axes._axes.Axes), "No or incorrect plot was created for transverse."
    plt.clf()

    results.breakthrough(x_position=1)
    assert isinstance(plt.gca(), matplotlib.axes._axes.Axes), "No or incorrect plot was created for breakthrough."
    plt.clf()

    results.plume_2d()
    assert isinstance(plt.gca(), matplotlib.axes._axes.Axes), "No or incorrect plot was created for plume_2d."
    plt.clf()

    results.plume_3d()
    assert isinstance(plt.gca(), matplotlib.axes._axes.Axes), "No or incorrect plot was created for plume_3d."
    plt.clf()
