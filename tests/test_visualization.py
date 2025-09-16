from types import NoneType
import pytest
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
import prettytable
from mibitrans.analysis.mass_balance import mass_balance
from mibitrans.visualize.plot_line import centerline
from mibitrans.visualize.plot_surface import plume_2d, plume_3d
from mibitrans.visualize.show_mass_balance import generate_mass_balance_tables
from mibitrans.transport.domenico import no_decay, linear_decay, instant_reaction
from tests.test_example_data import test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars

model_no_decay = no_decay(test_hydro_pars, test_ads_pars, test_source_pars, test_model_pars)
model_linear_decay = linear_decay(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars)
model_instant_reaction = instant_reaction(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars)
def test_centerline():
    centerline(model_no_decay)
    assert isinstance(plt.gca(), plt.Axes)

def test_plume_2d():
    plume_2d(model_no_decay)
    assert isinstance(plt.gca(), plt.Axes)

def test_plume_3d():
    ax = plume_3d(model_no_decay)
    assert isinstance(ax, mpl_toolkits.mplot3d.axes3d.Axes3D)

@pytest.mark.parametrize(
    "plottable, expected",
    [
        (-1, TypeError),
        (test_hydro_pars, TypeError),
    ]
)

def test_input_plotting(plottable, expected):
    with pytest.raises(expected):
        centerline(plottable)
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
        (3, TypeError)
    ])

def test_show_mass_balance(model, expected):
    if isinstance(expected, list):
        mb = mass_balance(model)
        table_list = list(generate_mass_balance_tables(mb))
        for i, table in enumerate(table_list):
            assert isinstance(table, expected[i])
    else:
        with pytest.raises(expected):
            mass_balance(model)