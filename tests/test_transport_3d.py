import warnings
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from mibitrans import AttenuationParameters
from mibitrans import ElectronAcceptors
from mibitrans import HydrologicalParameters
from mibitrans import ModelParameters
from mibitrans import SourceParameters
from mibitrans.data.check_input import DomainValueError

# from mibitrans.data.parameters import AttenuationParameters
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
        # Extent of y-domain should be at least the size of total source zone
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
        ({"mass_ratios": [0.2, 0.4]}, None),
        ({"mass_ratios": [0.2, -0.3]}, DomainValueError),
        ({"mass_ratios": "nonsense"}, TypeError),
        ({"mass_ratios": [0.2]}, ValueError),
        ({"mass_ratios": [0.2, 0.3, 0.1]}, ValueError),
    ],
)
def test_chain_decay_setup(
    pars, expected, test_hydro_pars, test_att_pars_chain, test_source_pars_chain, test_model_pars
) -> None:
    """Test if chain decay method accepts correct data and raises the intended errors if not."""
    model_object = Transport3DConcrete(test_hydro_pars, test_att_pars_chain, test_source_pars_chain, test_model_pars)
    if expected is None:
        model_object.chain_decay(**pars)
        assert model_object.mode == "chain_decay", (
            f"Model mode should be 'chain_decay', but is '{model_object.mode}' instead."
        )
    elif expected in [DomainValueError, ValueError, TypeError]:
        with pytest.raises(expected):
            model_object.chain_decay(**pars)


def test_chain_decay_setup_dataclasses(
    test_hydro_pars, test_att_pars, test_att_pars_chain, test_source_pars, test_source_pars_chain, test_model_pars
) -> None:
    """Test if chain decay method raises ValueError if att and/or source pars do not contain chain decay information."""
    model_object = Transport3DConcrete(test_hydro_pars, test_att_pars, test_source_pars_chain, test_model_pars)
    with pytest.raises(ValueError):
        model_object.chain_decay(mass_ratios=[0.8, 0.7])

    model_object = Transport3DConcrete(test_hydro_pars, test_att_pars_chain, test_source_pars, test_model_pars)
    with pytest.raises(ValueError):
        model_object.chain_decay(mass_ratios=[0.8, 0.7])

    model_object = Transport3DConcrete(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    with pytest.raises(ValueError):
        model_object.chain_decay(mass_ratios=[0.8, 0.7])


def test_mode_switch(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars, test_electron_acceptor_pars):
    """Test if model correctly switches modes when using the mode property."""
    model_object = Transport3DConcrete(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    with pytest.raises(ValueError):
        # Model should raise error if attempting to switch to core-fringe model without providing parameters.
        model_object.mode = "core-fringe"
    model_object.mode = "linear"
    assert model_object.mode == "linear", "Model mode was not switched to 'linear'."
    with pytest.raises(ValueError):
        # Model should raise error if attempting to switch to instant reaction for any other model than Bioscreen
        model_object.mode = "instant_reaction"


@pytest.mark.parametrize(
    "hydro, att, src, model, ea, expected",
    [
        (
            HydrologicalParameters(velocity=1, porosity=0.3, alpha_x=1, alpha_y=0.2, alpha_z=0.001),
            AttenuationParameters(decay_rate=0.001),
            SourceParameters(source_zone_boundary=5, source_zone_concentration=1, depth=3, total_mass=np.inf),
            ModelParameters(model_length=200, model_width=20, model_time=20 * 365, dx=1, dy=1, dt=365),
            ElectronAcceptors(0.5, 1, 1),
            [None, None],
        ),
        (
            HydrologicalParameters(velocity=1, porosity=0.3, alpha_x=1, alpha_y=0.2, alpha_z=0.001),
            AttenuationParameters(decay_rate=0.001),
            SourceParameters(source_zone_boundary=5, source_zone_concentration=1, depth=3, total_mass=np.inf),
            ModelParameters(model_length=2000, model_width=20, model_time=20 * 365, dx=1, dy=1, dt=365),
            ElectronAcceptors(0.5, 1, 1),
            [UserWarning, "Model domain is oversized for the estimated plume length"],
        ),
        (
            HydrologicalParameters(velocity=1, porosity=0.3, alpha_x=1, alpha_y=0.2, alpha_z=0.001),
            AttenuationParameters(decay_rate=0.001),
            SourceParameters(source_zone_boundary=5, source_zone_concentration=1, depth=3, total_mass=np.inf),
            ModelParameters(model_length=10, model_width=20, model_time=20 * 365, dx=1, dy=1, dt=365),
            ElectronAcceptors(0.5, 1, 1),
            [UserWarning, "Model domain is undersized for the estimated plume length"],
        ),
        (
            HydrologicalParameters(velocity=1, porosity=0.3, alpha_x=1, alpha_y=0.2, alpha_z=0.001),
            AttenuationParameters(decay_rate=0),
            SourceParameters(source_zone_boundary=5, source_zone_concentration=1, depth=3, total_mass=np.inf),
            ModelParameters(model_length=200, model_width=20, model_time=20 * 365, dx=1, dy=1, dt=365),
            None,
            [UserWarning, "Model domain is undersized for the estimated location of advective front"],
        ),
        (
            HydrologicalParameters(velocity=1, porosity=0.3, alpha_x=1, alpha_y=0.2, alpha_z=0.001),
            AttenuationParameters(decay_rate=0.001),
            SourceParameters(source_zone_boundary=5, source_zone_concentration=1, depth=3, total_mass=np.inf),
            ModelParameters(model_length=2000, model_width=20, model_time=365, dx=1, dy=1, dt=365),
            None,
            [UserWarning, "Model domain is oversized for the estimated location of advective front"],
        ),
    ],
)
def test_domain_check(hydro, att, src, model, ea, expected):
    """Test if domain check gives the correct warning if model domain is too small or too long."""
    model_obj = Transport3DConcrete(hydro, att, src, model)
    if ea is not None:
        model_obj.fringe_degradation(ea, 1)
    if expected[0] is UserWarning:
        with pytest.warns(UserWarning, match=expected[1]):
            model_obj.domain_check()
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            model_obj.domain_check()


def test_model_results_independance(test_mibitrans_model_fringe):
    """Test to make sure that changing model parameters does not change result parameters."""
    model, results = test_mibitrans_model_fringe
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
