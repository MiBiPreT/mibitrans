import copy
import numpy as np
import pytest
from mibitrans.data.parameters import AttenuationParameters
from mibitrans.data.parameters import HydrologicalParameters
from mibitrans.data.parameters import ModelParameters
from mibitrans.data.parameters import SourceParameters
from mibitrans.transport.models import Anatrans
from mibitrans.transport.models import Bioscreen
from mibitrans.transport.models import Mibitrans


@pytest.fixture(scope="session")
def test_hydro_pars():
    """HydrologicalParameters fixture with example data for tests."""
    return HydrologicalParameters(velocity=10 / 365, porosity=0.25, alpha_x=10, alpha_y=1, alpha_z=0.1)


@pytest.fixture(scope="session")
def test_att_pars():
    """AttenuationParameters fixture with example data for tests."""
    return AttenuationParameters(
        retardation=1,
        half_life=0.1 * 365,
    )


@pytest.fixture(scope="session")
def test_att_pars_nodecay():
    """AttenuationParameters fixture with example data for tests."""
    return AttenuationParameters(
        retardation=1,
        decay_rate=0,
    )


electron_acceptor_dict = dict(delta_oxygen=0.5, delta_nitrate=0.5, ferrous_iron=0.5, delta_sulfate=0.5, methane=0.5)


@pytest.fixture(scope="session")
def test_source_pars():
    """SourceParameters fixture with example data for tests."""
    return SourceParameters(
        source_zone_boundary=np.array([5, 10, 15]),
        source_zone_concentration=np.array([10, 5, 2]),
        depth=10,
        total_mass=1000000,
    )


@pytest.fixture(scope="session")
def test_model_pars():
    """ModelParameters fixture with example data for tests."""
    return ModelParameters(model_length=50, model_width=30, model_time=3 * 365, dx=10, dy=5, dt=1 * 365)


@pytest.fixture(scope="session")
def test_model_pars_short(test_source_pars, test_model_pars):
    """Model Parameters fixture with smaller model width for testing."""
    short_model_pars = copy.copy(test_model_pars)
    short_model_pars.model_width = test_source_pars.source_zone_boundary[-1] - 1
    return short_model_pars


@pytest.fixture(scope="session")
def test_mibitrans_model_nodecay(test_hydro_pars, test_att_pars_nodecay, test_source_pars, test_model_pars):
    """Mibitrans fixture model object for testing, with no decay."""
    obj = Mibitrans(test_hydro_pars, test_att_pars_nodecay, test_source_pars, test_model_pars)
    res = obj.run()
    return obj, res


@pytest.fixture(scope="session")
def test_mibitrans_model_lineardecay(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """Mibitrans fixture model object for testing, with linear decay."""
    obj = Mibitrans(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    res = obj.run()
    return obj, res


@pytest.fixture(scope="session")
def test_mibitrans_model_instantreaction(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """Mibitrans fixture model object for testing, with instant reaction."""
    obj = Mibitrans(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    obj.instant_reaction(electron_acceptors=electron_acceptor_dict)
    res = obj.run()
    return obj, res


@pytest.fixture(scope="module")
def test_anatrans_model_nodecay(test_hydro_pars, test_att_pars_nodecay, test_source_pars, test_model_pars):
    """Anatrans fixture model object for testing, with no decay."""
    obj = Anatrans(test_hydro_pars, test_att_pars_nodecay, test_source_pars, test_model_pars)
    res = obj.run()
    return obj, res


@pytest.fixture(scope="module")
def test_anatrans_model_lineardecay(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """Anatrans fixture model object for testing, with linear decay."""
    obj = Anatrans(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    res = obj.run()
    return obj, res


@pytest.fixture(scope="module")
def test_anatrans_model_instantreaction(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """Anatrans fixture model object for testing, with instant reaction."""
    obj = Anatrans(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    obj.instant_reaction(electron_acceptors=electron_acceptor_dict)
    res = obj.run()
    return obj, res


@pytest.fixture(scope="session")
def test_bioscreen_model_nodecay(test_hydro_pars, test_att_pars_nodecay, test_source_pars, test_model_pars):
    """Bioscreen fixture model object for testing, with no decay."""
    test_att_pars_nodecay.decay_rate = 0
    obj = Bioscreen(test_hydro_pars, test_att_pars_nodecay, test_source_pars, test_model_pars)
    res = obj.run()
    return obj, res


@pytest.fixture(scope="session")
def test_bioscreen_model_lineardecay(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """Bioscreen fixture model object for testing, with linear decay."""
    obj = Bioscreen(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    res = obj.run()
    return obj, res


@pytest.fixture(scope="session")
def test_bioscreen_model_instantreaction(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars):
    """Bioscreen fixture model object for testing, with instant reaction."""
    obj = Bioscreen(test_hydro_pars, test_att_pars, test_source_pars, test_model_pars)
    obj.instant_reaction(electron_acceptors=electron_acceptor_dict)
    res = obj.run()
    return obj, res
