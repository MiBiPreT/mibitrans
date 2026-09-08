"""Author: Jorrit Bakker.

Module handling testing of data input functionality
"""

import matplotlib.pyplot as plt
import numpy as np
import pytest
from mibitrans.data.check_input import DomainValueError
from mibitrans.data.check_input import MissingValueError
from mibitrans.data.parameter_information import FringeElectronAcceptors
from mibitrans.data.parameter_information import UtilizationFactor
from mibitrans.data.parameters import AttenuationParameters
from mibitrans.data.parameters import HydrologicalParameters
from mibitrans.data.parameters import ModelParameters
from mibitrans.data.parameters import SourceParameters


# Test HydrologicalParameters
@pytest.mark.parametrize(
    "parameters, error",
    [
        (dict(velocity=1, porosity=0.2, alpha_x=1, alpha_y=1), None),
        (dict(h_gradient=1, h_conductivity=1, porosity=0.2, alpha_x=1, alpha_y=1, alpha_z=0.1, diffusion=0.0004), None),
        (dict(), MissingValueError),
        (dict(porosity=0.2, alpha_x=1, alpha_y=1), MissingValueError),
        (dict(velocity=1, alpha_x=1, alpha_y=1), MissingValueError),
        (dict(velocity=1, porosity=0.2, alpha_y=1), MissingValueError),
        (dict(velocity=1, porosity=0.2, alpha_x=1), MissingValueError),
        (dict(h_gradient=1, porosity=0.2, alpha_x=1, alpha_y=1), MissingValueError),
        (dict(h_conductivity=1, porosity=0.2, alpha_x=1, alpha_y=1), MissingValueError),
        (dict(velocity="1", porosity=0.2, alpha_x=1, alpha_y=1), TypeError),
        (dict(velocity=1, porosity="2", alpha_x=1, alpha_y=1), TypeError),
        (dict(velocity=-1, porosity=0.2, alpha_x=1, alpha_y=1), DomainValueError),
        (dict(velocity=1, porosity=2, alpha_x=1, alpha_y=1), DomainValueError),
        (dict(velocity=1, porosity=0.2, alpha_x=1, alpha_y=1, alpha_z=-1), DomainValueError),
        (dict(velocity=1, porosity=0.2, alpha_x=1, alpha_y=1, diffusion=-0.1), DomainValueError),
    ],
)
def test_hyrologicalparameters_validation(parameters, error) -> None:
    """Test validation check of HydrologicalParameters dataclass."""
    if error is None:
        HydrologicalParameters(**parameters)
    else:
        with pytest.raises(error):
            HydrologicalParameters(**parameters)


@pytest.mark.parametrize(
    "parameters, value, error",
    [
        (dict(velocity=1, porosity=0.2, alpha_x=1, alpha_y=1), 2, None),
        (dict(velocity=1, porosity=0.2, alpha_x=1, alpha_y=1), "no", TypeError),
        (dict(velocity=1, porosity=0.2, alpha_x=1, alpha_y=1), -1, DomainValueError),
    ],
)
def test_hydrologicalparameters_setattribute(parameters, value, error) -> None:
    """Test validation of parameters after initialization."""
    hydro = HydrologicalParameters(**parameters)
    if error is None:
        hydro.alpha_x = value
    else:
        with pytest.raises(error):
            hydro.alpha_x = value


@pytest.mark.parametrize(
    "test, param, expected",
    [
        (dict(velocity=1, porosity=0.5, alpha_x=2, alpha_y=3), "velocity", 1),
        (dict(velocity=1, porosity=0.5, h_conductivity=2, h_gradient=2, alpha_x=2, alpha_y=3), "velocity", 8),
    ],
)
def test_hydrologicalparameters_output(test, param, expected) -> None:
    """Test output of HydrologicalParameters dataclass."""
    if "velocity" in test.keys() and "h_gradient" in test.keys():
        with pytest.warns(UserWarning):
            hydro = HydrologicalParameters(**test)
    else:
        hydro = HydrologicalParameters(**test)
    assert getattr(hydro, param) == expected


# Test AttenuationParameters
@pytest.mark.parametrize(
    "parameters, error",
    [
        (dict(decay_rate=0.2), None),
        (dict(half_life=2), None),
        (dict(retardation=1), None),
        (dict(decay_rate=0.2, retardation=2), None),
        (dict(bulk_density=1, partition_coefficient=1, fraction_organic_carbon=1), None),
        (dict(bulk_density=1, partition_coefficient=1, fraction_organic_carbon=1), None),
        (dict(), None),
        (dict(decay_rate=[0.2, 0.1, 0.05]), None),
        (dict(half_life=np.array([100, 73, 11])), None),
        (dict(decay_rate=[0.3]), None),
        (dict(decay_rate=1, half_life=1), UserWarning),
        (dict(half_life="one"), TypeError),
        (dict(half_life=-1), DomainValueError),
        (dict(retardation="one"), TypeError),
        (dict(retardation=0.1), DomainValueError),
        (dict(retardation=1, fraction_organic_carbon="no"), TypeError),
        (dict(retardation=1, fraction_organic_carbon=2), DomainValueError),
        (dict(decay_rate=[-2, 3]), DomainValueError),
        (dict(half_life=[300, [200, 100]]), TypeError),
    ],
)
def test_attenuationparameters_validation(parameters, error) -> None:
    """Test validation check of AdsorptionParameters dataclass."""
    if error is None:
        AttenuationParameters(**parameters)
    elif error is UserWarning:
        with pytest.warns(error):
            AttenuationParameters(**parameters)
    else:
        with pytest.raises(error):
            AttenuationParameters(**parameters)


@pytest.mark.parametrize(
    "test, param, expected",
    [
        (dict(retardation=1), "retardation", 1),
        (dict(decay_rate=2), "decay_rate", 2),
        (dict(decay_rate=2), "half_life", np.log(2) / 2),
        (dict(half_life=2), "half_life", 2),
        (dict(half_life=2), "decay_rate", np.log(2) / 2),
        (dict(decay_rate=2, half_life=2), "decay_rate", 2),
        (dict(decay_rate=2, half_life=2), "half_life", np.log(2) / 2),
        (dict(decay_rate=[2, 4]), "half_life", [np.log(2) / 2, np.log(2) / 4]),
    ],
)
@pytest.mark.filterwarnings("ignore:Both contaminant decay rate")
def test_attenuationparameters_output(test, param, expected) -> None:
    """Test output of AttenuationParameters dataclass."""
    att = AttenuationParameters(**test)
    assert getattr(att, param) == pytest.approx(expected)


@pytest.mark.parametrize(
    "test, parameter, value, error",
    [
        (dict(retardation=1), "retardation", 1, None),
        (
            dict(bulk_density=1, partition_coefficient=1, fraction_organic_carbon=1),
            "fraction_organic_carbon",
            {},
            TypeError,
        ),
        (
            dict(bulk_density=1, partition_coefficient=1, fraction_organic_carbon=1),
            "fraction_organic_carbon",
            2,
            DomainValueError,
        ),
        (dict(half_life=3), "half_life", 0, None),
    ],
)
def test_attenuationparameters_setattribute(test, value, parameter, error) -> None:
    """Test validation of parameters after initialization."""
    att = AttenuationParameters(**test)
    if error is None:
        setattr(att, parameter, value)
        assert getattr(att, parameter) == value
    else:
        with pytest.raises(error):
            setattr(att, parameter, value)


@pytest.mark.parametrize(
    "test, expected",
    [
        (dict(util_oxygen=1, util_nitrate=2, util_ferrous_iron=3, util_sulfate=4, util_methane=5), None),
        (dict(util_oxygen="1", util_nitrate=2, util_ferrous_iron=3, util_sulfate=4, util_methane=5), TypeError),
        (dict(util_oxygen=-1, util_nitrate=2, util_ferrous_iron=3, util_sulfate=4, util_methane=5), DomainValueError),
    ],
)
def test_attenuationparameters_utilization(test, expected, test_att_pars) -> None:
    """Test set_utilization_factor method of AttenuationParameters dataclass."""
    if expected is None:
        UtilizationFactor(**test)
    else:
        with pytest.raises(expected):
            UtilizationFactor(**test)


# Test SourceParameters
@pytest.mark.parametrize(
    "parameters, error",
    [
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5), None),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5, total_mass=2), None),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5, total_mass="inf"), None),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5, total_mass="infint"), None),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5, total_mass=np.inf), None),
        (dict(source_zone_boundary=1, source_zone_concentration=[3], depth=5, total_mass=2), None),
        (
            dict(source_zone_boundary=np.array([1, 2, 3]), source_zone_concentration=[3, 2, 1], depth=5, total_mass=2),
            None,
        ),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[[3, 2], [4, 3], [6, 5]], depth=5), None),
        (dict(source_zone_boundary=5, source_zone_concentration=[5, 4, 3, 2, 1], depth=5), None),
        (dict(source_zone_boundary=(1, 2), source_zone_concentration=[3, 2], depth=5, total_mass=2), TypeError),
        (dict(source_zone_boundary=["one", 2], source_zone_concentration=[3, 2], depth=5, total_mass=2), TypeError),
        (
            dict(source_zone_boundary=[1, 2], source_zone_concentration=np.array([3, -2]), depth=5, total_mass=2),
            DomainValueError,
        ),
        (
            dict(source_zone_boundary=np.array(([1, 2, 3], [4, 5, 6])), source_zone_concentration=[3, 2], depth=5),
            ValueError,
        ),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[2, 3], depth=5, total_mass=2), ValueError),
        (dict(source_zone_boundary=[-1, 2], source_zone_concentration=[3, 2], depth=5, total_mass=2), DomainValueError),
        (dict(source_zone_boundary=-1, source_zone_concentration=3), DomainValueError),
        (dict(source_zone_boundary=[1, 2, 3], source_zone_concentration=[3, 2], depth=5, total_mass=2), ValueError),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth="five", total_mass=2), TypeError),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=-5, total_mass=2), DomainValueError),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5, total_mass=[2, 3]), TypeError),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5, total_mass=-2), DomainValueError),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5, total_mass="nons"), ValueError),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[[3, 2], [4, 3], 1], depth=5), TypeError),
        (
            dict(source_zone_boundary=[1, 2, 3], source_zone_concentration=np.array([[3, 2, 1], [4, 3, 2]]), depth=5),
            ValueError,
        ),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[[3, 2, 1], [4, 3, 2]], depth=5), ValueError),
        (
            dict(source_zone_boundary=[1, 2, 3], source_zone_concentration=[[3, 2, 1], [4, 3, -2]], depth=5),
            DomainValueError,
        ),
    ],
)
def test_sourceparameters_validation(parameters, error) -> None:
    """Test validation check of SourceParameters dataclass."""
    if error is None:
        SourceParameters(**parameters)
    else:
        with pytest.raises(error):
            SourceParameters(**parameters)


# Test SourceParameters
@pytest.mark.parametrize(
    "parameter, value, error",
    [
        ("source_zone_boundary", [1, 2, 3], None),
        ("source_zone_concentration", [3, 2, 1], None),
        ("total_mass", 1000, None),
        ("total_mass", "infini", None),
        ("total_mass", np.inf, None),
        ("source_zone_concentration", [[3, 2, 1], [4, 3, 2]], None),
        ("source_zone_concentration", [1, 2, 3], ValueError),
        ("source_zone_concentration", 1, ValueError),
        ("source_zone_concentration", "No", TypeError),
        ("source_zone_concentration", [[3, 2], [4, 3]], ValueError),
    ],
)
def test_sourceparameters_validation_setattr(parameter, value, error) -> None:
    """Test validation check of SourceParameters dataclass. when setting attributes."""
    src = SourceParameters(
        source_zone_boundary=np.array([1, 2, 3]), source_zone_concentration=np.array([3, 2, 1]), depth=5
    )
    if error is None:
        setattr(src, parameter, value)
    else:
        with pytest.raises(error):
            setattr(src, parameter, value)


def test_sourceparameters_visualize_zone():
    """Test if source zone visualization creates a plot."""
    source = SourceParameters(np.array([1, 2, 3]), np.array([3, 2, 1]), 10, 1000)
    source.visualize_source_zone()
    assert isinstance(plt.gca(), plt.Axes)
    plt.clf()
    source = SourceParameters(np.array([1, 2, 3]), [np.array([6, 3, 1]), np.array([3, 2, 1])], 10)
    source.visualize_source_zone()
    assert isinstance(plt.gca(), plt.Axes)


def test_sourceparameters_visualize_depletion():
    """Test if source depletion visualization creates a plot."""
    source = SourceParameters(np.array([1, 2, 3]), np.array([3, 2, 1]), 10, 1000)
    hydro = HydrologicalParameters(velocity=1, porosity=0.2, alpha_x=1, alpha_y=1)
    source.visualize_source_depletion(hydro)
    assert isinstance(plt.gca(), plt.Axes)
    source = SourceParameters(np.array([1, 2, 3]), [np.array([3, 2, 1]), np.array([4, 3, 2])], 10, 1000)
    with pytest.raises(Exception):
        source.visualize_source_depletion(hydro)
    source = SourceParameters(np.array([1, 2, 3]), [np.array([3, 2, 1]), np.array([4, 3, 2])], 10, np.inf)
    with pytest.raises(ValueError):
        source.visualize_source_depletion(hydro)


# Test ModelParameters
@pytest.mark.parametrize(
    "parameters, error",
    [
        (dict(model_length=1, model_width=1, model_time=1), None),
        (dict(model_length=1, model_width=1, model_time=1, dx=1, dy=1, dt=1), None),
        (dict(dx=1, dy=1, dt=1), TypeError),
        (dict(model_length="one", model_width=1, model_time=1), TypeError),
        (dict(model_length=-2, model_width=1, model_time=1), DomainValueError),
        (dict(model_length=1, model_width=1, model_time=1, dx=2), ValueError),
        (dict(model_length=1, model_width=1, model_time=1, dy=2), ValueError),
        (dict(model_length=1, model_width=1, model_time=1, dt=2), ValueError),
    ],
)
def test_modelparameters_validation(parameters, error) -> None:
    """Test validation check of ModelParameters dataclass."""
    if error is None:
        ModelParameters(**parameters)
    else:
        with pytest.raises(error):
            ModelParameters(**parameters)


@pytest.mark.parametrize(
    "parameter, value, error",
    [
        ("model_length", 3, None),
        ("dy", 0.2, None),
        ("model_time", 3, None),
        ("model_length", "nonsense", TypeError),
        ("model_length", 0.1, ValueError),
        ("dx", 2, ValueError),
        ("model_width", "nonsense", TypeError),
        ("model_width", 0.1, ValueError),
        ("dy", 2, ValueError),
        ("model_time", "nonsense", TypeError),
        ("model_time", 0.1, ValueError),
        ("dt", 2, ValueError),
    ],
)
def test_modelparameters_validation_setattr(parameter, value, error) -> None:
    """Test validation check of ModelParameters dataclass when setting attributes."""
    mod = ModelParameters(model_length=1, model_width=1, model_time=1, dx=0.5, dy=0.5, dt=0.5)
    if error is None:
        setattr(mod, parameter, value)
    else:
        with pytest.raises(error):
            setattr(mod, parameter, value)


@pytest.mark.parametrize(
    "test, param, expected",
    [
        (dict(model_length=1, model_width=1, model_time=1, dx=0.5), "model_length", 1),
        (dict(model_length=1, model_width=1, model_time=1, dx=0.5), "dx", 0.5),
    ],
)
def test_modelparameters_output(test, param, expected) -> None:
    """Test output of ModelParameters dataclass."""
    model = ModelParameters(**test)
    assert model.__dict__[param] == expected


def test_calculation_optional_discretization():
    """Test if model discretization is calculated if not given."""
    model = ModelParameters(model_length=100, model_width=50, model_time=10)
    assert model.dx, "Model dx should have been calculated, but was not."
    assert model.dy, "Model dy should have been calculated, but was not."
    assert model.dt, "Model dt should have been calculated, but was not."


@pytest.mark.parametrize(
    "parameters, error",
    [
        (
            dict(electron_acceptor_concentration=5, stoichiometric_ratio=3.4, electron_acceptor_molecular_weight=54),
            None,
        ),
        (
            dict(
                electron_acceptor_concentration=[5, 4, 3],
                stoichiometric_ratio=[3.4, 5.6, 2.1],
                electron_acceptor_molecular_weight=[34.2, 67, 52.1],
            ),
            None,
        ),
        (
            dict(
                electron_acceptor_concentration=[5, 4, 3],
                stoichiometric_ratio=[3.4, 5.6, 2.1],
                electron_acceptor_molecular_weight=54,
            ),
            ValueError,
        ),
        (
            dict(
                electron_acceptor_concentration=[5, 4, 3],
                stoichiometric_ratio=[3.4, -5.6, 2.1],
                electron_acceptor_molecular_weight=[34.2, 67, 52.1],
            ),
            DomainValueError,
        ),
        (
            dict(
                electron_acceptor_concentration="der",
                stoichiometric_ratio=[3.4, 5.6, 2.1],
                electron_acceptor_molecular_weight=[34.2, 67, 52.1],
            ),
            TypeError,
        ),
        (
            dict(
                electron_acceptor_concentration=[5, "die", 3],
                stoichiometric_ratio=[3.4, 5.6, 2.1],
                electron_acceptor_molecular_weight=[34.2, 67, 52.1],
            ),
            TypeError,
        ),
    ],
)
def test_fringe_electron_acceptor(parameters, error) -> None:
    """Test validation of FringeElectronAcceptor dataclass."""
    if error is None:
        FringeElectronAcceptors(**parameters)
    else:
        with pytest.raises(error):
            FringeElectronAcceptors(**parameters)


@pytest.mark.parametrize(
    "parameters, expected",
    [
        (dict(ea=FringeElectronAcceptors(6, 1, 1), electron_donor_molecular_weight=1), 6),
        (dict(ea=FringeElectronAcceptors(5, 7.5, 32), electron_donor_molecular_weight=78.11), 1.6272916666666666),
        (
            dict(ea=FringeElectronAcceptors([5, 6], [7.5, 6], [32, 62]), electron_donor_molecular_weight=78.11),
            2.8871303763440856,
        ),
        (dict(ea=FringeElectronAcceptors(6, 1, 1), electron_donor_molecular_weight="heavy"), TypeError),
        (dict(ea=FringeElectronAcceptors(6, 1, 1), electron_donor_molecular_weight=[5, 4, 3]), TypeError),
        (dict(ea=FringeElectronAcceptors(6, 1, 1), electron_donor_molecular_weight=-10), DomainValueError),
    ],
)
def test_fringe_electron_acceptor_calculate_bc(parameters, expected) -> None:
    """Test calculation of biodegradation capacity from electron acceptors."""
    ea = parameters["ea"]
    if isinstance(expected, (int, float)):
        bc = ea.calculate_bc(parameters["electron_donor_molecular_weight"])
        assert bc == expected
    else:
        with pytest.raises(expected):
            ea.calculate_bc(parameters["electron_donor_molecular_weight"])
