import numpy as np
import pytest
from mibitrans.transport.domenico import Domenico, no_decay, linear_decay, instant_reaction

from tests.test_example_data import testingdata_nodecay, testingdata_lineardecay, testingdata_instantreaction
from tests.test_example_data import test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars

@pytest.mark.parametrize(
    "hydro, ads, source, model",
    [
        (test_hydro_pars, test_ads_pars, test_source_pars, test_model_pars),
    ])

def test_domenico_parent(hydro, ads, source, model) -> None:
    parent = Domenico(hydro, ads, source, model)
    shape_arrays = (len(parent.t), len(parent.y), len(parent.x))
    # Source zone concentrations adapted for superposition should still have the same length as those in input
    assert (len(parent.c_source) == len(source.source_zone_concentration)) and (len(parent.c_source) == len(source.source_zone_boundary))
    # Extent of y-domain should be at least the size of
    assert (np.max(parent.y) + abs(np.min(parent.y))) >= (np.max(source.source_zone_boundary) * 2)
    assert parent.xxx.shape == shape_arrays
    assert parent.yyy.shape == shape_arrays
    assert parent.ttt.shape == shape_arrays
    assert parent.ads_pars.retardation is not None
    assert hydro.velocity / parent.ads_pars.retardation == parent.rv

@pytest.mark.parametrize(
    "model, expected",
    [
        (no_decay(test_hydro_pars, test_ads_pars, test_source_pars, test_model_pars), testingdata_nodecay),
        (linear_decay(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars), testingdata_lineardecay),
        (instant_reaction(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars), testingdata_instantreaction)
    ])

def test_transport_equations_numerical(model, expected):
    assert model.cxyt == pytest.approx(expected)
