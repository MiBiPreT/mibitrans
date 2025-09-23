"""Author: Jorrit Bakker.

Module handling testing of data input functionality
"""

import numpy as np
import pytest
from mibitrans.data.read import AdsorptionParameters
from mibitrans.data.read import DegradationParameters
from mibitrans.data.read import HydrologicalParameters
from mibitrans.data.read import ModelParameters
from mibitrans.data.read import SourceParameters


# Test HydrologicalParameters
@pytest.mark.parametrize(
    "parameters, error",
    [
        (dict(velocity=1, porosity=0.2, alpha_x=1, alpha_y=1), None),
        (dict(h_gradient=1, h_conductivity=1, porosity=0.2, alpha_x=1, alpha_y=1), None),
        (dict(), ValueError),
        (dict(porosity=0.2, alpha_x=1, alpha_y=1), ValueError),
        (dict(velocity=1, alpha_x=1, alpha_y=1), ValueError),
        (dict(velocity=1, porosity=0.2, alpha_y=1), ValueError),
        (dict(velocity=1, porosity=0.2, alpha_x=1), ValueError),
        (dict(h_gradient=1, porosity=0.2, alpha_x=1, alpha_y=1), ValueError),
        (dict(h_conductivity=1, porosity=0.2, alpha_x=1, alpha_y=1), ValueError),
        (dict(velocity="1", porosity=0.2, alpha_x=1, alpha_y=1), TypeError),
        (dict(velocity=1, porosity="2", alpha_x=1, alpha_y=1), TypeError),
        (dict(velocity=-1, porosity=0.2, alpha_x=1, alpha_y=1), ValueError),
        (dict(velocity=1, porosity=2, alpha_x=1, alpha_y=1), ValueError),
        (dict(velocity=1, porosity=0.2, alpha_x=1, alpha_y=1, alpha_z=-1), ValueError),
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
        (dict(velocity=1, porosity=0.2, alpha_x=1, alpha_y=1), -1, ValueError),
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
def test_hyrologicalparameters_output(test, param, expected) -> None:
    """Test output of HydrologicalParameters dataclass."""
    if "velocity" in test.keys() and "h_gradient" in test.keys():
        with pytest.warns(UserWarning):
            hydro = HydrologicalParameters(**test)
    else:
        hydro = HydrologicalParameters(**test)
    assert getattr(hydro, param) == expected


# Test AdsorptionParameters
@pytest.mark.parametrize(
    "parameters, error",
    [
        (dict(retardation=1), None),
        (dict(bulk_density=1, partition_coefficient=1, fraction_organic_carbon=1), None),
        (dict(), ValueError),
        (dict(bulk_density=1, partition_coefficient=1), ValueError),
        (dict(retardation="one"), TypeError),
        (dict(retardation=0.1), ValueError),
        (dict(retardation=1, fraction_organic_carbon="no"), TypeError),
        (dict(retardation=1, fraction_organic_carbon=2), ValueError),
    ],
)
def test_adsorptionparameters_validation(parameters, error) -> None:
    """Test validation check of AdsorptionParameters dataclass."""
    if error is None:
        AdsorptionParameters(**parameters)
    else:
        with pytest.raises(error):
            AdsorptionParameters(**parameters)


@pytest.mark.parametrize(
    "parameters, value, error",
    [
        (dict(bulk_density=1, partition_coefficient=1, fraction_organic_carbon=1), 0.01, None),
        (dict(bulk_density=1, partition_coefficient=1, fraction_organic_carbon=1), {}, TypeError),
        (dict(bulk_density=1, partition_coefficient=1, fraction_organic_carbon=1), 2, ValueError),
    ],
)
def test_adsorptionparameters_setattribute(parameters, value, error) -> None:
    """Test validation of parameters after initialization."""
    ads = AdsorptionParameters(**parameters)
    if error is None:
        ads.fraction_organic_carbon = value
    else:
        with pytest.raises(error):
            ads.fraction_organic_carbon = value


@pytest.mark.parametrize(
    "test, param, expected",
    [
        (dict(retardation=1), "retardation", 1),
    ],
)
def test_adsorptionparameters_output(test, param, expected) -> None:
    """Test output of AdsorptionParameters dataclass."""
    ads = AdsorptionParameters(**test)
    assert getattr(ads, param) == expected


# Test DegradationParameters
@pytest.mark.parametrize(
    "parameters, error",
    [
        (dict(decay_rate=0.2), None),
        (dict(delta_oxygen=1, delta_nitrate=1, ferrous_iron=1, delta_sulfate=1, methane=1), None),
        (dict(), ValueError),
        (dict(delta_oxygen=1, delta_nitrate=1, ferrous_iron=1, delta_sulfate=1), ValueError),
        (dict(half_life="one"), TypeError),
        (dict(half_life=-1), ValueError),
    ],
)
def test_degradationparameters_validation(parameters, error) -> None:
    """Test validation check of DegradationParameters dataclass."""
    if error is None:
        DegradationParameters(**parameters)
    else:
        with pytest.raises(error):
            DegradationParameters(**parameters)


@pytest.mark.parametrize(
    "parameters, value, error",
    [
        (dict(delta_oxygen=1, delta_nitrate=1, ferrous_iron=1, delta_sulfate=1, methane=1), 2, None),
        (dict(delta_oxygen=1, delta_nitrate=1, ferrous_iron=1, delta_sulfate=1, methane=1), "dont", TypeError),
        (dict(delta_oxygen=1, delta_nitrate=1, ferrous_iron=1, delta_sulfate=1, methane=1), -1, ValueError),
    ],
)
def test_degradationparameters_setattribute(parameters, value, error) -> None:
    """Test validation of parameters after initialization."""
    deg = DegradationParameters(**parameters)
    if error is None:
        deg.delta_oxygen = value
    else:
        with pytest.raises(error):
            deg.delta_oxygen = value


@pytest.mark.parametrize(
    "test, param, expected",
    [
        (dict(half_life=2), "decay_rate", np.log(2) / 2),
        (dict(decay_rate=0.5, half_life=2), "decay_rate", 0.5),
    ],
)
def test_degradationparameters_output(test, param, expected) -> None:
    """Test output of DegradationParameters dataclass."""
    if "half_life" in test.keys() and "decay_rate" in test.keys():
        with pytest.warns(UserWarning):
            deg = DegradationParameters(**test)
    else:
        deg = DegradationParameters(**test)
    assert getattr(deg, param) == expected


deg_test_object = DegradationParameters(delta_oxygen=1, delta_nitrate=1, ferrous_iron=1, delta_sulfate=1, methane=1)


@pytest.mark.parametrize(
    "test, expected",
    [
        (dict(util_oxygen=1, util_nitrate=2, util_ferrous_iron=3, util_sulfate=4, util_methane=5), None),
        (dict(util_oxygen="1", util_nitrate=2, util_ferrous_iron=3, util_sulfate=4, util_methane=5), TypeError),
        (dict(util_oxygen=-1, util_nitrate=2, util_ferrous_iron=3, util_sulfate=4, util_methane=5), ValueError),
    ],
)
def test_degradationparameters_utilization(test, expected) -> None:
    """Test set_utilization_factor method of DegradationParameters dataclass."""
    if expected is None:
        deg_test_object.set_utilization_factor(**test)
    else:
        with pytest.raises(expected):
            deg_test_object.set_utilization_factor(**test)


def test_degradationparameters_utilization_setattr() -> None:
    """Test attribute setting of utilization_factor of DegradationParameters dataclass."""
    deg_test_object.utilization_factor = deg_test_object.utilization_factor
    with pytest.raises(TypeError):
        deg_test_object.utilization_factor = "not a dataclass"


# Test SourceParameters
@pytest.mark.parametrize(
    "parameters, error",
    [
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5), None),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5, total_mass=2), None),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5, total_mass="inf"), None),
        (dict(source_zone_boundary=1, source_zone_concentration=[3], depth=5, total_mass=2), None),
        (
            dict(source_zone_boundary=np.array([1, 2, 3]), source_zone_concentration=[3, 2, 1], depth=5, total_mass=2),
            None,
        ),
        (dict(), ValueError),
        (dict(source_zone_boundary=(1, 2), source_zone_concentration=[3, 2], depth=5, total_mass=2), TypeError),
        (dict(source_zone_boundary=["one", 2], source_zone_concentration=[3, 2], depth=5, total_mass=2), TypeError),
        (
            dict(source_zone_boundary=[1, 2], source_zone_concentration=np.array([-3, 2]), depth=5, total_mass=2),
            ValueError,
        ),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[2, 3], depth=5, total_mass=2), ValueError),
        (dict(source_zone_boundary=[-1, 2], source_zone_concentration=[3, 2], depth=5, total_mass=2), ValueError),
        (dict(source_zone_boundary=-1, source_zone_concentration=3), ValueError),
        (dict(source_zone_boundary=[1, 2, 3], source_zone_concentration=[3, 2], depth=5, total_mass=2), ValueError),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth="five", total_mass=2), TypeError),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=-5, total_mass=2), ValueError),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5, total_mass=[2, 3]), TypeError),
        (dict(source_zone_boundary=[1, 2], source_zone_concentration=[3, 2], depth=5, total_mass=-2), ValueError),
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
        ("source_zone_concentration", [1, 2, 3], ValueError),
        ("source_zone_concentration", 1, ValueError),
        ("source_zone_concentration", "No", TypeError),
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


@pytest.mark.parametrize(
    "test, param, expected",
    [
        (
            dict(source_zone_boundary=[1, 2, 3], source_zone_concentration=[6, 4, 2], depth=5, total_mass=2),
            "source_zone_boundary",
            np.array([1, 2, 3]),
        ),
        (
            dict(source_zone_boundary=[2, 3, 1], source_zone_concentration=[4, 2, 6], depth=5, total_mass=2),
            "source_zone_boundary",
            np.array([1, 2, 3]),
        ),
        (
            dict(source_zone_boundary=[2, 3, 1], source_zone_concentration=[4, 2, 6], depth=5, total_mass=2),
            "source_zone_concentration",
            np.array([6, 4, 2]),
        ),
        (
            dict(source_zone_boundary=[1, 2, 3], source_zone_concentration=[6, 4, 2], depth=5, total_mass="inf"),
            "total_mass",
            "infinite",
        ),
    ],
)
def test_sourceparameters_output(test, param, expected) -> None:
    """Test output of SourceParameters dataclass."""
    unordered = np.array(test["source_zone_boundary"]) < test["source_zone_boundary"][0]
    if True in unordered:
        with pytest.warns(UserWarning):
            source = SourceParameters(**test)
    else:
        source = SourceParameters(**test)
    assert source.__dict__[param] == pytest.approx(expected)


# Test ModelParameters
@pytest.mark.parametrize(
    "parameters, error",
    [
        (dict(), None),
        (dict(model_length=1, model_width=1, model_time=1, dx=1, dy=1, dt=1), None),
        (dict(model_length="one"), TypeError),
        (dict(model_length=-2), ValueError),
        (dict(model_length=1, dx=2), ValueError),
        (dict(model_width=1, dy=2), ValueError),
        (dict(model_time=1, dt=2), ValueError),
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
        (dict(model_length=1, dx=0.5), "model_length", 1),
        (dict(model_length=1, dx=0.5), "dx", 0.5),
    ],
)
def test_modelparameters_output(test, param, expected) -> None:
    """Test output of ModelParameters dataclass."""
    model = ModelParameters(**test)
    assert model.__dict__[param] == expected
