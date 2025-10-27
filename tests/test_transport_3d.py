import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from mibitrans.data.parameters import AttenuationParameters
from mibitrans.transport.model_parent import Transport3D


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
        parent = Transport3D(*args)
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
            parent = Transport3D(*args)
            range_y = abs(parent.y[0]) + abs(parent.y[-1])
            assert range_y >= parent._src_pars.source_zone_boundary[-1] * 2
    elif error is TypeError:
        with pytest.raises(error):
            Transport3D(*args)


@pytest.mark.parametrize(
    "att, expected",
    [
        (AttenuationParameters(retardation=1), 1),
        (AttenuationParameters(bulk_density=2, partition_coefficient=10, fraction_organic_carbon=0.03), 3.4),
    ],
)
def test_retardation_calculation(att, expected, test_hydro_pars, test_source_pars, test_model_pars) -> None:
    """Test if retardation is calculated correctly when Transport3D class is initialized."""
    parent = Transport3D(test_hydro_pars, att, test_source_pars, test_model_pars)
    assert parent._att_pars.retardation == expected


def test_plotting_methods(test_domenico_nodecay_model):
    """Test if plotting methods defined in parent model are working."""
    test_domenico_nodecay_model.centerline()
    assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
    plt.clf()

    test_domenico_nodecay_model.transverse(x_position=1)
    assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
    plt.clf()

    test_domenico_nodecay_model.breakthrough(x_position=1)
    assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
    plt.clf()

    test_domenico_nodecay_model.plume_2d()
    assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
    plt.clf()

    test_domenico_nodecay_model.plume_3d()
    assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
    plt.clf()
