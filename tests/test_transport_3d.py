import numpy as np
import pytest
from mibitrans.data.parameters import AttenuationParameters
from mibitrans.data.parameters import ModelParameters
from mibitrans.transport.model_parent import Transport3D
from tests.test_example_data import test_att_pars
from tests.test_example_data import test_deg_pars
from tests.test_example_data import test_hydro_pars
from tests.test_example_data import test_model_pars
from tests.test_example_data import test_source_pars

short_width_model_pars = ModelParameters(model_length=50, model_width=30, model_time=3 * 365, dx=10, dy=5, dt=1 * 365)
short_width_model_pars.model_width = test_source_pars.source_zone_boundary[-1] - 1


@pytest.mark.parametrize(
    "hydro, att, source, model, error",
    [
        (test_hydro_pars, test_att_pars, test_source_pars, test_model_pars, None),
        (1, test_att_pars, test_source_pars, test_model_pars, TypeError),
        (test_hydro_pars, test_hydro_pars, test_source_pars, test_model_pars, TypeError),
        (test_hydro_pars, test_att_pars, "wrong", test_model_pars, TypeError),
        (test_hydro_pars, test_att_pars, test_source_pars, test_deg_pars, TypeError),
        (test_hydro_pars, test_att_pars, test_source_pars, short_width_model_pars, UserWarning),
    ],
)
@pytest.mark.filterwarnings("ignore:UserWarning")
def test_transport_parent(hydro, att, source, model, error) -> None:
    """Test functionality, results and errors of Transport3D parent class."""
    if error is None:
        parent = Transport3D(hydro, att, source, model)
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
        assert parent.att_pars.retardation is not None
        assert hydro.velocity / parent.att_pars.retardation == parent.rv
    elif error is UserWarning:
        with pytest.warns(UserWarning):
            parent = Transport3D(hydro, att, source, model)
            range_y = abs(parent.y[0]) + abs(parent.y[-1])
            assert range_y >= parent.src_pars.source_zone_boundary[-1] * 2
    elif error is TypeError:
        with pytest.raises(error):
            Transport3D(hydro, att, source, model)


@pytest.mark.parametrize(
    "att, expected",
    [
        (AttenuationParameters(retardation=1), 1),
        (AttenuationParameters(bulk_density=2, partition_coefficient=10, fraction_organic_carbon=0.03), 3.4),
    ],
)
def test_retardation_calculation(att, expected):
    """Test if retardation is calculated correctly when Transport3D class is initialized."""
    parent = Transport3D(test_hydro_pars, att, test_source_pars, test_model_pars)
    assert parent.att_pars.retardation == expected
