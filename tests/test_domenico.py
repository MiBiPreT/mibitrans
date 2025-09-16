import numpy as np
import pytest
from mibitrans.transport.domenico import Domenico
from mibitrans.transport.domenico import instant_reaction
from mibitrans.transport.domenico import linear_decay
from mibitrans.transport.domenico import no_decay
from tests.test_example_data import test_ads_pars
from tests.test_example_data import test_deg_pars
from tests.test_example_data import test_hydro_pars
from tests.test_example_data import test_model_pars
from tests.test_example_data import test_source_pars
from tests.test_example_data import testingdata_instantreaction
from tests.test_example_data import testingdata_lineardecay
from tests.test_example_data import testingdata_nodecay


@pytest.mark.parametrize(
    "hydro, ads, source, model, error",
    [
        (test_hydro_pars, test_ads_pars, test_source_pars, test_model_pars, None),
        (1, test_ads_pars, test_source_pars, test_model_pars, TypeError),
        (test_hydro_pars, test_hydro_pars, test_source_pars, test_model_pars, TypeError),
        (test_hydro_pars, test_ads_pars, "wrong", test_model_pars, TypeError),
        (test_hydro_pars, test_ads_pars, test_source_pars, test_deg_pars, TypeError),
    ],
)
def test_domenico_parent(hydro, ads, source, model, error) -> None:
    """Test functionality, results and errors of Domenico parent class."""
    if error is None:
        parent = Domenico(hydro, ads, source, model)
        shape_arrays = (len(parent.t), len(parent.y), len(parent.x))
        # Source zone concentrations adapted for superposition should still have the same length as those in input
        assert (len(parent.c_source) == len(source.source_zone_concentration)) and (
            len(parent.c_source) == len(source.source_zone_boundary)
        )
        # Extent of y-domain should be at least the size of
        assert (np.max(parent.y) + abs(np.min(parent.y))) >= (np.max(source.source_zone_boundary) * 2)
        assert parent.xxx.shape == shape_arrays
        assert parent.yyy.shape == shape_arrays
        assert parent.ttt.shape == shape_arrays
        assert parent.ads_pars.retardation is not None
        assert hydro.velocity / parent.ads_pars.retardation == parent.rv
    else:
        with pytest.raises(error):
            parent = Domenico(hydro, ads, source, model)


@pytest.mark.parametrize(
    "model, expected",
    [
        (no_decay(test_hydro_pars, test_ads_pars, test_source_pars, test_model_pars), testingdata_nodecay),
        (
            linear_decay(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars),
            testingdata_lineardecay,
        ),
        (
            instant_reaction(test_hydro_pars, test_ads_pars, test_deg_pars, test_source_pars, test_model_pars),
            testingdata_instantreaction,
        ),
    ],
)
def test_transport_equations_numerical(model, expected):
    """Test numerical output of transport equation child classes of Domenico."""
    assert model.cxyt == pytest.approx(expected)
