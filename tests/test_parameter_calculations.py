"""Author: Jorrit Bakker.

File testing functionality of parameter_calculations module.
"""
import numpy as np
import pytest
from anatrans.analysis.parameter_calculations import calculate_biodegradation_capacity
from anatrans.analysis.parameter_calculations import calculate_dispersivity
from anatrans.analysis.parameter_calculations import calculate_flow_velocity
from anatrans.analysis.parameter_calculations import calculate_linear_decay
from anatrans.analysis.parameter_calculations import calculate_retardation
from anatrans.analysis.parameter_calculations import calculate_source_decay
from anatrans.analysis.parameter_calculations import calculate_source_decay_instant


@pytest.mark.parametrize(
    "test, expected",
    [
        (dict(R=1), 1),
        (dict(R=1, rho=2, n=0.3, Koc=4, foc=5),  1),
        (dict(rho=2, n=0.5, Koc=20, foc=1e-5),  1.0008),
    ])

def test_calculate_retardation(test, expected):
    """Test calculation of retardation factor."""
    assert calculate_retardation(test) == pytest.approx(expected)

@pytest.mark.parametrize(
    "test, expected",
    [
        (dict(v=1), 1),
        (dict(v=1, k=10, i=0.005, n = 0.5),  1),
        (dict(k=10, i=0.005, n = 0.5),  0.1),
    ])

def test_calculate_flow_velocity(test, expected):
    """Test calculation of flow velocity."""
    assert calculate_flow_velocity(test) == pytest.approx(expected)

@pytest.mark.parametrize(
    "test, expected",
    [
        (dict(alpha_x=1, alpha_y=1, alpha_z=1), (1, 1, 1)),
        (dict(alpha_x=1, alpha_y=1, alpha_z=1, lp=1),  (1, 1, 1)),
        (dict(lp=127.092354820933),  (5, 0.5, 0.05)),
        (dict(alpha_x=1, alpha_y=1, alpha_z=0),  (1, 1, 1e-10)),
    ])

def test_calculate_dispersivity(test, expected):
    """Test calculation of dispersivity."""
    assert calculate_dispersivity(test) == pytest.approx(expected)

@pytest.mark.parametrize(
    "test, expected",
    [
        (dict(mu=0.1), 0.1),
        (dict(mu=0.1, t_half=5),  0.1),
        (dict(t_half=5),  np.log(2) / 5),
    ])

def test_calculate_linear_decay(test, expected):
    """Test calculation of linear decay coefficient."""
    assert calculate_linear_decay(test) == pytest.approx(expected)

@pytest.mark.parametrize(
    "test, expected",
    [
        (dict(c_source=np.array([[0,10], [10,5], [30,2], [50,0]]), v=1, n=0.5, d_source=10, m_total=100), 0.024),
    ])

def test_calculate_source_decay(test, expected):
    """Test calculation of source decay coefficient."""
    assert calculate_source_decay(test) == pytest.approx(expected)

@pytest.mark.parametrize(
    "test, expected",
    [
        (dict(dO=1.65, dNO3=0.7, Fe2=16.6, dSO4=22.4, CH4=6.6, c_source=np.array([[0,10], [10,5], [30,2], [50,0]]),
              v=1, n=0.5, d_source=10, m_total=100), 0.0972864932405937),
    ])

def test_calculate_source_decay_instant(test, expected):
    """Test calculation of source decay coefficient for instant reaction model."""
    biodeg_cap = calculate_biodegradation_capacity(test)
    assert calculate_source_decay_instant(test, biodeg_cap) == pytest.approx(expected)

@pytest.mark.parametrize(
    "test, expected",
    [
        (dict(dO=1.65, dNO3=0.7, Fe2=16.6, dSO4=22.4, CH4=6.6), 14.657298648118742),
    ])

def test_calculate_biodegradation_capacity(test, expected):
    """Test calculation of biodegradation capacity."""
    assert calculate_biodegradation_capacity(test) == pytest.approx(expected)
