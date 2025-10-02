from types import NoneType
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
import numpy as np
import prettytable
import pytest
from mibitrans.analysis.mass_balance import mass_balance
from mibitrans.data.read import SourceParameters
from mibitrans.transport.domenico import InstantReaction
from mibitrans.transport.domenico import LinearDecay
from mibitrans.transport.domenico import NoDecay
from mibitrans.visualize.plot_line import breakthrough
from mibitrans.visualize.plot_line import centerline
from mibitrans.visualize.plot_line import transverse
from mibitrans.visualize.plot_surface import plume_2d
from mibitrans.visualize.plot_surface import plume_3d
from mibitrans.visualize.show_conditions import source_zone
from mibitrans.visualize.show_mass_balance import generate_mass_balance_tables
from tests.test_example_data import test_ads_pars
from tests.test_example_data import test_deg_pars
from tests.test_example_data import test_hydro_pars
from tests.test_example_data import test_model_pars
from tests.test_example_data import test_source_pars

model_no_decay = NoDecay(test_hydro_pars, test_ads_pars, test_source_pars, test_model_pars)
model_linear_decay = LinearDecay(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars)
model_instant_reaction = InstantReaction(
    test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars
)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, plt.Axes),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
def test_centerline(animate, expected):
    """Test if plot object is generated in centerline function."""
    if not animate:
        centerline(model_no_decay, animate=animate)
        assert isinstance(plt.gca(), expected)
    else:
        ani = centerline(model_no_decay, animate=animate)
        assert isinstance(ani, expected)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, plt.Axes),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
def test_transverse(animate, expected):
    """Test if plot object is generated in transverse function."""
    if not animate:
        transverse(model_no_decay, x_position=1, animate=animate)
        assert isinstance(plt.gca(), expected)
    else:
        ani = transverse(model_no_decay, x_position=1, animate=animate)
        assert isinstance(ani, expected)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, plt.Axes),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
def test_breakthrough(animate, expected):
    """Test if plot object is generated in breakthrough function."""
    if not animate:
        breakthrough(model_no_decay, x_position=1, animate=animate)
        assert isinstance(plt.gca(), plt.Axes)
    else:
        ani = breakthrough(model_no_decay, x_position=1, animate=animate)
        assert isinstance(ani, expected)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, plt.Axes),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
def test_plume_2d(animate, expected):
    """Test if plot object is generated in plume 2d function."""
    if not animate:
        plume_2d(model_no_decay, animate=animate)
        assert isinstance(plt.gca(), expected)
    else:
        ani = plume_2d(model_no_decay, animate=animate)
        assert isinstance(ani, expected)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, mpl_toolkits.mplot3d.axes3d.Axes3D),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
def test_plume_3d(animate, expected):
    """Test if plot object is generated in plume 3d function."""
    ax = plume_3d(model_no_decay, animate=animate)
    assert isinstance(ax, expected)


source = SourceParameters(np.array([1, 2, 3]), np.array([3, 2, 1]), 10, 1000)


def test_source_zone():
    """Test if plot object is generated in source zone function."""
    source_zone(source)
    assert isinstance(plt.gca(), plt.Axes)


@pytest.mark.parametrize(
    "plottable, expected",
    [
        (-1, TypeError),
        (test_hydro_pars, TypeError),
    ],
)
def test_input_plotting(plottable, expected):
    """Test if input validation plotting function."""
    with pytest.raises(expected):
        centerline(plottable)
    with pytest.raises(expected):
        transverse(plottable, x_position=1)
    with pytest.raises(expected):
        breakthrough(plottable, x_position=1)
    with pytest.raises(expected):
        plume_2d(plottable)
    with pytest.raises(expected):
        plume_3d(plottable)


pt_type = prettytable.prettytable.PrettyTable


@pytest.mark.parametrize(
    "model, expected",
    [
        (model_no_decay, [pt_type, pt_type, NoneType, NoneType, NoneType]),
        (model_linear_decay, [pt_type, pt_type, pt_type, NoneType, NoneType]),
        (model_instant_reaction, [pt_type, pt_type, NoneType, pt_type, pt_type]),
        (3, TypeError),
    ],
)
def test_show_mass_balance(model, expected):
    """Test show_mass_balance function input validation and output generation."""
    if isinstance(expected, list):
        mb = mass_balance(model)
        table_list = list(generate_mass_balance_tables(mb))
        for i, table in enumerate(table_list):
            assert isinstance(table, expected[i])
    else:
        with pytest.raises(expected):
            generate_mass_balance_tables(model)


@pytest.mark.parametrize(
    "model",
    [
        model_linear_decay,
        model_instant_reaction,
    ],
)
def test_print_mass_balance(model):
    """Test if print_mass_balance runs without error."""
    mb = mass_balance(model)
    generate_mass_balance_tables(mb)
