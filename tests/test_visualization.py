import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
import pytest
from mibitrans.data.check_input import DomainValueError
from mibitrans.transport.models import Mibitrans
from mibitrans.visualize.plot_line import breakthrough
from mibitrans.visualize.plot_line import centerline
from mibitrans.visualize.plot_line import transverse
from mibitrans.visualize.plot_surface import plume_2d
from mibitrans.visualize.plot_surface import plume_3d

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
def test_centerline(animate, expected, test_anatrans_model_nodecay, test_anatrans_model_lineardecay):
    """Test if plot object is generated in centerline function."""
    model_nd, results_nd = test_anatrans_model_nodecay
    model_lin, results_lin = test_anatrans_model_lineardecay
    if not animate:
        centerline(results_nd, animate=animate)
        assert isinstance(plt.gca(), expected)
        plt.clf()
        centerline(
            [results_nd, results_lin],
            legend_names=["no decay", "linear decay"],
            animate=animate,
        )
        assert isinstance(plt.gca(), expected)

    else:
        ani = centerline(results_nd, animate=animate)
        assert isinstance(ani, expected)
        anim = centerline(
            [results_nd, results_lin],
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
def test_parameter_check_centerline(y_pos, time, expected, test_anatrans_model_nodecay):
    """Test if centerline function properly raises warnings and errors for parameters outside domain."""
    model, results = test_anatrans_model_nodecay
    if isinstance(expected, matplotlib.axes._axes.Axes):
        centerline(results, y_pos, time)
        assert isinstance(plt.gca(), expected)
    elif expected is UserWarning:
        with pytest.warns(expected):
            centerline(results, y_pos, time)
            assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
    elif expected is TypeError or expected is DomainValueError:
        with pytest.raises(expected):
            centerline(results, y_pos, time)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, matplotlib.axes._axes.Axes),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
@pytest.mark.filterwarnings("ignore:Animation was deleted")
def test_transverse(animate, expected, test_anatrans_model_nodecay, test_anatrans_model_lineardecay):
    """Test if plot object is generated in transverse function."""
    model_nd, results_nd = test_anatrans_model_nodecay
    model_lin, results_lin = test_anatrans_model_lineardecay
    if not animate:
        transverse(results_nd, x_position=10, animate=animate)
        assert isinstance(plt.gca(), expected)
        plt.clf()
        transverse(
            [results_nd, results_lin],
            x_position=10,
            legend_names=["no decay", "linear decay"],
            animate=animate,
        )
        assert isinstance(plt.gca(), expected)

    else:
        ani = transverse(results_nd, x_position=10, animate=animate)
        assert isinstance(ani, expected)
        anim = transverse(
            [results_nd, results_lin],
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
def test_parameter_check_transverse(x_pos, time, expected, test_anatrans_model_nodecay):
    """Test if transverse function properly raises warnings and errors for parameters outside domain."""
    model, results = test_anatrans_model_nodecay
    if isinstance(expected, matplotlib.axes._axes.Axes):
        transverse(results, x_pos, time)
        assert isinstance(plt.gca(), expected)
    elif expected is UserWarning:
        with pytest.warns(expected):
            transverse(results, x_pos, time)
            assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
    elif expected is TypeError or expected is DomainValueError:
        with pytest.raises(expected):
            transverse(results, x_pos, time)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, matplotlib.axes._axes.Axes),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
@pytest.mark.filterwarnings("ignore:Animation was deleted")
def test_breakthrough(animate, expected, test_anatrans_model_nodecay, test_anatrans_model_lineardecay):
    """Test if plot object is generated in breakthrough function."""
    model_nd, results_nd = test_anatrans_model_nodecay
    model_lin, results_lin = test_anatrans_model_lineardecay
    if not animate:
        breakthrough(results_nd, x_position=10, animate=animate)
        assert isinstance(plt.gca(), expected)
        plt.clf()
        breakthrough(
            [results_nd, results_lin],
            x_position=10,
            legend_names=["no decay", "linear decay"],
            animate=animate,
        )
        assert isinstance(plt.gca(), expected)
    else:
        ani = breakthrough(results_nd, x_position=10, animate=animate)
        assert isinstance(ani, expected)
        anim = breakthrough(
            [results_nd, results_lin],
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
def test_parameter_check_breakthrough(x_pos, y_pos, expected, test_anatrans_model_nodecay):
    """Test if breakthrough function properly raises warnings and errors for parameters outside domain."""
    model, results = test_anatrans_model_nodecay
    if isinstance(expected, matplotlib.axes._axes.Axes):
        breakthrough(results, x_pos, y_pos)
        assert isinstance(plt.gca(), expected)
    elif expected is UserWarning:
        with pytest.warns(expected):
            breakthrough(results, x_pos, y_pos)
            assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
    elif expected is TypeError or expected is DomainValueError:
        with pytest.raises(expected):
            breakthrough(results, x_pos, y_pos)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, matplotlib.axes._axes.Axes),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
@pytest.mark.filterwarnings("ignore:Animation was deleted")
def test_plume_2d(animate, expected, test_anatrans_model_nodecay):
    """Test if plot object is generated in plume 2d function."""
    model, results = test_anatrans_model_nodecay
    if not animate:
        plume_2d(results, animate=animate)
        assert isinstance(plt.gca(), expected)
    else:
        ani = plume_2d(results, animate=animate)
        assert isinstance(ani, expected)


@pytest.mark.parametrize(
    "animate, expected",
    [
        (False, mpl_toolkits.mplot3d.axes3d.Axes3D),
        (True, matplotlib.animation.FuncAnimation),
    ],
)
@pytest.mark.filterwarnings("ignore:Animation was deleted")
def test_plume_3d(animate, expected, test_anatrans_model_nodecay):
    """Test if plot object is generated in plume 3d function."""
    model, results = test_anatrans_model_nodecay
    ax = plume_3d(results, animate=animate)
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


@pytest.mark.parametrize(
    "plotter, expected",
    [
        ("centerline", matplotlib.axes._axes.Axes),
        ("transverse", matplotlib.axes._axes.Axes),
        ("breakthrough", matplotlib.axes._axes.Axes),
        ("plume_2d", matplotlib.axes._axes.Axes),
        ("plume_3d", matplotlib.axes._axes.Axes),
        ("plume_3d", mpl_toolkits.mplot3d.axes3d.Axes3D),
    ],
)
def test_innate_plotting_methods(plotter, expected, test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """Test if the plotting class methods of the results object work."""
    model = Mibitrans(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    results = model.run()
    match plotter:
        case "centerline":
            results.centerline()
        case "transverse":
            results.transverse(1)
        case "breakthrough":
            results.breakthrough(1)
        case "plume_2d":
            results.plume_2d()
        case "plume_3d":
            ax = results.plume_3d()
            assert isinstance(ax, expected)
    assert isinstance(plt.gca(), expected)
    plt.clf()
