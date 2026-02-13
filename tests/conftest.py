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

# Test parameters loosely based on Keesler site. Some adaptations to allow for more robust tests.


@pytest.fixture(scope="session")
def test_hydro_pars():
    """HydrologicalParameters fixture with example data for tests."""
    return HydrologicalParameters(
        h_gradient=0.048,  # [m/m]
        h_conductivity=0.495,  # [m/d]
        porosity=0.25,  # [-]
        alpha_x=4.1,  # [m]
        alpha_y=0.4,  # [m]
        alpha_z=0.01,  # [m]
    )


@pytest.fixture(scope="session")
def test_att_pars():
    """AttenuationParameters fixture with example data for tests."""
    return AttenuationParameters(
        bulk_density=1.7,  # [kg/L]
        partition_coefficient=38,  # [L/kg]
        fraction_organic_carbon=0.000057,  # [-]
        half_life=0.5 * 365,  # [1/day]
    )


@pytest.fixture(scope="session")
def test_att_pars_nodecay():
    """AttenuationParameters fixture with example data for tests."""
    return AttenuationParameters(
        bulk_density=1.7,  # [kg/L]
        partition_coefficient=38,  # [L/kg]
        fraction_organic_carbon=0.000057,  # [-]
        decay_rate=0,  # [-]
    )


electron_acceptor_dict = dict(
    delta_oxygen=2.05 - 0.4,  # [g/m3]
    delta_nitrate=0.07 - 0,  # [g/m3]
    ferrous_iron=16.6,  # [g/m3]
    delta_sulfate=26.2 - 3.8,  # [g/m3]
    methane=6.6,  # [g/m3]
)


@pytest.fixture(scope="session")
def test_source_pars():
    """SourceParameters fixture with example data for tests."""
    return SourceParameters(
        source_zone_boundary=np.array([2, 11, 20]),  # [m]
        source_zone_concentration=np.array([13.68, 2.508, 0.057]),  # [g/m3]
        depth=3,  # [m]
        total_mass=2000000,  # [g]
    )


@pytest.fixture(scope="session")
def test_model_pars():
    """ModelParameters fixture with example data for tests."""
    return ModelParameters(
        model_length=100,  # [m]
        model_width=40,  # [m]
        model_time=5 * 365,  # [days]
        dx=20,  # [m]
        dy=10,  # [m]
        dt=365,  # [days]
    )


@pytest.fixture(scope="session")
def test_model_pars_short(test_source_pars, test_model_pars):
    """Model Parameters fixture with smaller model width for testing."""
    short_model_pars = copy.copy(test_model_pars)
    short_model_pars.model_width = test_source_pars.source_zone_boundary[-1]
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
