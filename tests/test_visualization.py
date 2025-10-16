from types import NoneType
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
import prettytable
import pytest
from mibitrans.analysis.mass_balance import mass_balance
from mibitrans.data.check_input import DomainValueError
from mibitrans.visualize.plot_line import breakthrough
from mibitrans.visualize.plot_line import centerline
from mibitrans.visualize.plot_line import transverse
from mibitrans.visualize.plot_surface import plume_2d
from mibitrans.visualize.plot_surface import plume_3d
from mibitrans.visualize.show_mass_balance import generate_mass_balance_tables

matplotlib.use("Agg")  # Fixes tkinter.TclError in local tests


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, matplotlib.axes._axes.Axes),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
# Ignore warning for animation being deleted without being shown, as behaviour is intentional for testing purposes.
@pytest.mark.filterwarnings("ignore:Animation was deleted")
def test_centerline(animate, expected, test_domenico_nodecay_model, test_domenico_lineardecay_model):
    """Test if plot object is generated in centerline function."""
    if not animate:
        centerline(test_domenico_nodecay_model, animate=animate)
        assert isinstance(plt.gca(), expected)
        plt.clf()
        centerline(
            [test_domenico_nodecay_model, test_domenico_lineardecay_model],
            legend_names=["no decay", "linear decay"],
            animate=animate,
        )
        assert isinstance(plt.gca(), expected)

    else:
        ani = centerline(test_domenico_nodecay_model, animate=animate)
        assert isinstance(ani, expected)
        anim = centerline(
            [test_domenico_nodecay_model, test_domenico_lineardecay_model],
            legend_names=["no decay", "linear decay"],
            animate=animate,
        )
        assert isinstance(anim, expected)


@pytest.mark.parametrize(
    "y_pos, time, expected",
    [
        (1, 365, matplotlib.axes._axes.Axes),
        (None, 365, matplotlib.axes._axes.Axes),
        (-1, None, matplotlib.axes._axes.Axes),
        (100000, 365, UserWarning),
        (1, 40000000, UserWarning),
        ("nonsense", 365, TypeError),
        (1, -1, DomainValueError),
    ],
)
def test_parameter_check_centerline(y_pos, time, expected, test_domenico_nodecay_model):
    """Test if centerline function properly raises warnings and errors for parameters outside domain."""
    if isinstance(expected, matplotlib.axes._axes.Axes):
        centerline(test_domenico_nodecay_model, y_pos, time)
        assert isinstance(plt.gca(), expected)
    elif expected is UserWarning:
        with pytest.warns(expected):
            centerline(test_domenico_nodecay_model, y_pos, time)
            assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
    elif expected is TypeError or expected is DomainValueError:
        with pytest.raises(expected):
            centerline(test_domenico_nodecay_model, y_pos, time)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, matplotlib.axes._axes.Axes),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
@pytest.mark.filterwarnings("ignore:Animation was deleted")
def test_transverse(animate, expected, test_domenico_nodecay_model, test_domenico_lineardecay_model):
    """Test if plot object is generated in transverse function."""
    if not animate:
        transverse(test_domenico_nodecay_model, x_position=10, animate=animate)
        assert isinstance(plt.gca(), expected)
        plt.clf()
        transverse(
            [test_domenico_nodecay_model, test_domenico_lineardecay_model],
            x_position=10,
            legend_names=["no decay", "linear decay"],
            animate=animate,
        )
        assert isinstance(plt.gca(), expected)

    else:
        ani = transverse(test_domenico_nodecay_model, x_position=10, animate=animate)
        assert isinstance(ani, expected)
        anim = transverse(
            [test_domenico_nodecay_model, test_domenico_lineardecay_model],
            x_position=10,
            legend_names=["no decay", "linear decay"],
            animate=animate,
        )
        assert isinstance(anim, expected)


@pytest.mark.parametrize(
    "x_pos, time, expected",
    [
        (10, 365, matplotlib.axes._axes.Axes),
        (-1, None, matplotlib.axes._axes.Axes),
        (1000000, 365, UserWarning),
        (10, 3000000, UserWarning),
        ("nonsense", 365, TypeError),
        (1, -1, DomainValueError),
    ],
)
def test_parameter_check_transverse(x_pos, time, expected, test_domenico_nodecay_model):
    """Test if transverse function properly raises warnings and errors for parameters outside domain."""
    if isinstance(expected, matplotlib.axes._axes.Axes):
        transverse(test_domenico_nodecay_model, x_pos, time)
        assert isinstance(plt.gca(), expected)
    elif expected is UserWarning:
        with pytest.warns(expected):
            transverse(test_domenico_nodecay_model, x_pos, time)
            assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
    elif expected is TypeError or expected is DomainValueError:
        with pytest.raises(expected):
            transverse(test_domenico_nodecay_model, x_pos, time)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, matplotlib.axes._axes.Axes),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
@pytest.mark.filterwarnings("ignore:Animation was deleted")
def test_breakthrough(animate, expected, test_domenico_nodecay_model, test_domenico_lineardecay_model):
    """Test if plot object is generated in breakthrough function."""
    if not animate:
        breakthrough(test_domenico_nodecay_model, x_position=10, animate=animate)
        assert isinstance(plt.gca(), expected)
        plt.clf()
        breakthrough(
            [test_domenico_nodecay_model, test_domenico_lineardecay_model],
            x_position=10,
            legend_names=["no decay", "linear decay"],
            animate=animate,
        )
        assert isinstance(plt.gca(), expected)
    else:
        ani = breakthrough(test_domenico_nodecay_model, x_position=10, animate=animate)
        assert isinstance(ani, expected)
        anim = breakthrough(
            [test_domenico_nodecay_model, test_domenico_lineardecay_model],
            x_position=10,
            legend_names=["no decay", "linear decay"],
            animate=animate,
        )
        assert isinstance(anim, expected)


@pytest.mark.parametrize(
    "x_pos, y_pos, expected",
    [
        (10, 1, matplotlib.axes._axes.Axes),
        (-1, None, matplotlib.axes._axes.Axes),
        (10000000, 365, UserWarning),
        (10, 200000000, UserWarning),
        ("nonsense", 365, TypeError),
        (-10, 1, DomainValueError),
    ],
)
def test_parameter_check_breakthrough(x_pos, y_pos, expected, test_domenico_nodecay_model):
    """Test if breakthrough function properly raises warnings and errors for parameters outside domain."""
    if isinstance(expected, matplotlib.axes._axes.Axes):
        breakthrough(test_domenico_nodecay_model, x_pos, y_pos)
        assert isinstance(plt.gca(), expected)
    elif expected is UserWarning:
        with pytest.warns(expected):
            breakthrough(test_domenico_nodecay_model, x_pos, y_pos)
            assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
    elif expected is TypeError or expected is DomainValueError:
        with pytest.raises(expected):
            breakthrough(test_domenico_nodecay_model, x_pos, y_pos)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, matplotlib.axes._axes.Axes),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
@pytest.mark.filterwarnings("ignore:Animation was deleted")
def test_plume_2d(animate, expected, test_domenico_nodecay_model):
    """Test if plot object is generated in plume 2d function."""
    if not animate:
        plume_2d(test_domenico_nodecay_model, animate=animate)
        assert isinstance(plt.gca(), expected)
    else:
        ani = plume_2d(test_domenico_nodecay_model, animate=animate)
        assert isinstance(ani, expected)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, mpl_toolkits.mplot3d.axes3d.Axes3D),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
@pytest.mark.filterwarnings("ignore:Animation was deleted")
def test_plume_3d(animate, expected, test_domenico_nodecay_model):
    """Test if plot object is generated in plume 3d function."""
    ax = plume_3d(test_domenico_nodecay_model, animate=animate)
    assert isinstance(ax, expected)


def test_source_zone(test_source_pars):
    """Test if plot object is generated in source zone function."""
    test_source_pars.visualize()
    assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)


@pytest.mark.parametrize(
    "plottable, expected",
    [
        (-1, TypeError),
        ("test_hydro_pars", TypeError),
    ],
)
def test_input_plotting(plottable, expected, request):
    """Test if input validation plotting function."""
    if isinstance(plottable, str):
        plottable = request.getfixturevalue(plottable)
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
        ("test_domenico_nodecay_model", [pt_type, pt_type, NoneType, NoneType, NoneType]),
        ("test_domenico_lineardecay_model", [pt_type, pt_type, pt_type, NoneType, NoneType]),
        ("test_domenico_instantreaction_model", [pt_type, pt_type, NoneType, pt_type, pt_type]),
        (3, TypeError),
    ],
)
@pytest.mark.filterwarnings("ignore:Decay rate was set")
def test_show_mass_balance(model, expected, request):
    """Test show_mass_balance function input validation and output generation."""
    if isinstance(model, str):
        model = request.getfixturevalue(model)
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
        "test_domenico_lineardecay_model",
        "test_domenico_instantreaction_model",
    ],
)
@pytest.mark.filterwarnings("ignore:UserWarning")
def test_print_mass_balance(model, request):
    """Test if print_mass_balance runs without error."""
    model_object = request.getfixturevalue(model)
    mb = mass_balance(model_object)
    generate_mass_balance_tables(mb)
