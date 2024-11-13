"""Author: Jorrit Bakker.

File testing functionality of parameter_calculations module.
"""

import pytest
from anatrans.analysis.parameter_calculations import calculate_dispersivity
from anatrans.analysis.parameter_calculations import calculate_flow_velocity
from anatrans.analysis.parameter_calculations import calculate_linear_decay
from anatrans.analysis.parameter_calculations import calculate_retardation
from anatrans.analysis.parameter_calculations import calculate_source_decay

@pytest.mark.parametrize(
    "test, expected",
    [
        (dict(R=1), 1),
        (dict(R=1, rho=2, n=0.3, Koc=4, foc=5),  1),
        (dict(rho=1.7, n=0.3, Koc=38, foc=5.7e-5),  1),
    ])

def test_calculate_retardation(test, expected):
    """Test calculation of retardation factor."""
    r = calculate_retardation(test)

    assert r == pytest.approx(expected, 0.02)
