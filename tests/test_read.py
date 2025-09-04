"""Author: Jorrit Bakker.

Module handling testing of data input functionality
"""

import pytest
import numpy as np
from mibitrans.data.read import HydrologicalParameters, AdsorptionDegradationParameters, ModelParameters, SourceParameters
from mibitrans.data.parameter_information import key_dictionary as k_dict

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
    ])

def test_hyrologicalparameters_validation(parameters, error) -> None:
    if error is None:
        HydrologicalParameters(**parameters)
    else:
        with pytest.raises(error):
            HydrologicalParameters(**parameters)

@pytest.mark.parametrize(
    "test, param, expected",
    [
        (dict(velocity=1, porosity=0.5, alpha_x=2, alpha_y=3), "velocity", 1),
        (dict(velocity=1, porosity=0.5, h_conductivity=2, h_gradient=2, alpha_x=2, alpha_y=3), "velocity", 8),
    ])

def test_hyrologicalparameters_output(test, param, expected) -> None:
    hydro = HydrologicalParameters(**test)
    assert hydro.__dict__[param] == expected

# Test AdsorptionDegradationParameters
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
        (dict(retardation=1, half_life="half"), TypeError),
        (dict(retardation=1, half_life=-1), ValueError),
    ])

def test_adsorptiondegradationparameters_validation(parameters, error) -> None:
    if error is None:
        AdsorptionDegradationParameters(**parameters)
    else:
        with pytest.raises(error):
            AdsorptionDegradationParameters(**parameters)

@pytest.mark.parametrize(
    "test, param, expected",
    [
        (dict(retardation=1), "retardation", 1),
        (dict(retardation=1, half_life=2), "decay_rate", np.log(2) / 2),
        (dict(retardation=1, decay_rate=0.5, half_life=2), "decay_rate", 0.5),
    ])

def test_adsorptiondegradationparameters_output(test, param, expected) -> None:
    ads_deg = AdsorptionDegradationParameters(**test)
    assert ads_deg.__dict__[param] == expected

# Test ModelParameters
@pytest.mark.parametrize(
    "parameters, error",
    [
        (dict(), None),
        (dict(model_length=1, model_width=1, model_time=1, dx=1, dy=1, dt=1), None),
        (dict(model_length="one"), TypeError),
        (dict(model_length=-2), ValueError),
    ])

def test_modelparameters_validation(parameters, error) -> None:
    if error is None:
        ModelParameters(**parameters)
    else:
        with pytest.raises(error):
            ModelParameters(**parameters)

@pytest.mark.parametrize(
    "test, param, expected",
    [
        (dict(model_length=1, dx=0.5), "model_length", 1),
        (dict(model_length=1, dx=0.5), "dx", 0.5),
        (dict(model_length=1, dx=2), "dx", 1),
    ])

def test_modelparameters_output(test, param, expected) -> None:
    model = ModelParameters(**test)
    assert model.__dict__[param] == expected

# Test SourceParameters
@pytest.mark.parametrize(
    "parameters, error",
    [
        (dict(source_zone_boundary=[1,2], source_zone_concentration=[3,2], depth=5), None),
        (dict(source_zone_boundary=[1,2], source_zone_concentration=[3,2], depth=5, total_mass=2), None),
        (dict(source_zone_boundary=[1,2], source_zone_concentration=[3,2], depth=5, total_mass="inf"), None),
        (dict(source_zone_boundary=1, source_zone_concentration=[3], depth=5, total_mass=2), None),
        (dict(source_zone_boundary=np.array([1,2,3]), source_zone_concentration=[3,2,1], depth=5, total_mass=2), None),
        (dict(), ValueError),
        (dict(source_zone_boundary=(1,2), source_zone_concentration=[3,2], depth=5, total_mass=2), TypeError),
        (dict(source_zone_boundary=["one",2], source_zone_concentration=[3,2], depth=5, total_mass=2), TypeError),
        (dict(source_zone_boundary=[1,2], source_zone_concentration=np.array([-3,2]), depth=5, total_mass=2), ValueError),
        (dict(source_zone_boundary=[-1,2], source_zone_concentration=[3,2], depth=5, total_mass=2), ValueError),
        (dict(source_zone_boundary=-1, source_zone_concentration=3), ValueError),
        (dict(source_zone_boundary=[1, 2, 3], source_zone_concentration=[3, 2], depth=5, total_mass=2), ValueError),
        (dict(source_zone_boundary=[1,2], source_zone_concentration=[3,2], depth="five", total_mass=2), TypeError),
        (dict(source_zone_boundary=[1,2], source_zone_concentration=[3,2], depth=-5, total_mass=2), ValueError),
        (dict(source_zone_boundary=[1,2], source_zone_concentration=[3,2], depth=5, total_mass=[2,3]), TypeError),
        (dict(source_zone_boundary=[1,2], source_zone_concentration=[3,2], depth=5, total_mass=-2), ValueError),
    ])

def test_sourceparameters_validation(parameters, error) -> None:
    if error is None:
        SourceParameters(**parameters)
    else:
        with pytest.raises(error):
            SourceParameters(**parameters)

@pytest.mark.parametrize(
    "test, param, expected",
    [
        (dict(source_zone_boundary=[1,2,3], source_zone_concentration=[6,4,2], depth=5, total_mass=2), "source_zone_boundary", np.array([1,2,3])),
        (dict(source_zone_boundary=[2,3,1], source_zone_concentration=[4,2,6], depth=5, total_mass=2), "source_zone_boundary", np.array([1,2,3])),
        (dict(source_zone_boundary=[2,3,1], source_zone_concentration=[4,2,6], depth=5, total_mass=2), "source_zone_concentration", np.array([6,4,2])),
        (dict(source_zone_boundary=[1,2,3], source_zone_concentration=[6,4,2], depth=5, total_mass="inf"), "total_mass", "infinite"),
        (dict(source_zone_boundary=[1,2,3], source_zone_concentration=[6,4,2], depth=5, total_mass="nonsense"),"total_mass", "infinite"),
    ])

def test_sourceparameters_output(test, param, expected) -> None:
    source = SourceParameters(**test)
    assert source.__dict__[param] == pytest.approx(expected)

##################
##### Legacy #####
##################

@pytest.mark.parametrize(
    "test, param, expected",
    [
        (dict(model_length=1, dx=0.5), "model_length", 1),
        (dict(model_length=1, dx=0.5), "dx", 0.5),
        (dict(model_length=1, dx=2), "dx", 1),
    ])

def test_sourceparameters_output(test, param, expected) -> None:
    source = SourceParameters(**test)
    assert source.__dict__[param] == expected

@pytest.mark.parametrize(
    "test, expected",
    [
        ({k_dict["v"][0] : 1, k_dict["R"][0] : 2, k_dict["mu"][0] : 3}, dict(v = 1, R = 2, mu = 3)),
        (dict(v = 1, R = 2, nonsense = 3), dict(v = 1, R = 2)),
        (dict(nonsense = 3), dict()),
        (dict(), dict()),
        ({k_dict["v"][-1] : 1}, {"v" : 1})
    ])

def test_from_dict(test, expected) -> None:
    """Test if from_dict gives expected output for various input dictionaries."""
    result = rd.from_dict(test)
    assert result == expected
