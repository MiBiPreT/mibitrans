"""Author: Jorrit Bakker.

Module handling data input in the form of a dictionary.
"""

import numpy as np
import pytest
import anatrans.data.check_input as ci


@pytest.mark.parametrize(
    "test, mode, expected",
    [
        (dict(v=1, lp=2, R=3, d_source=4, c_source=5, m_total=6, n=0.5), None, True),
        (dict(v=1, lp=2, R=3, c_source=5, m_total=6, n=0.5), None, False),
        (dict(v=1, lp=2, R=3, d_source=4, m_total=6, n=0.5), None, False),
        (dict(v=1, lp=2, R=3, d_source=4, c_source=5, n=0.5), None, False),
        (dict(v=1, lp=2, R=3, d_source=4, c_source=5, m_total=6), None, False),
        (dict(), None, False),
        (dict(k=0.5, i=1, n=0.5, lp=2, R=3, d_source=4, c_source=5, m_total=6), None, True),
        (dict(k=0.5, n=1.5, lp=2, R=3, d_source=4, c_source=5, m_total=6), None, False),
        (dict(v=1, alpha_x=1.5, alpha_y=2, alpha_z=2.5, R=3, d_source=4, c_source=5, m_total=6, n=0.5), None, True),
        (dict(v=1, lp=2, rho=2.5, Koc=3, foc=3.5, n=0.5, d_source=4, c_source=5, m_total=6), None, True),
        (dict(k=0.5, i=1, n=0.5, lp=2, R=3, d_source=4, c_source=5, m_total=6, mu=7), "linear_decay", True),
        (dict(k=0.5, i=1, n=0.5, lp=2, R=3, d_source=4, c_source=5, m_total=6), "linear_decay", False),
        (dict(k=0.5, i=1, n=1.5, lp=2, R=3, d_source=4, c_source=5, m_total=6, dO=7, dNO3=8, Fe2=9, dSO4=10, CH4=11),
         "instant_reaction", True),
        (dict(k=0.5, i=1, n=0.5, lp=2, R=3, d_source=4, c_source=5, m_total=6), "instant_reaction", False),
        (dict(k=0.5, i=1, n=0.5, lp=2, R=3, d_source=5, m_total=6, dO=7, dNO3=8, Fe2=9, CH4=11),
         "instant_reaction", False),
    ])

def test_check_parameter(test, mode, expected) -> None:
    """Test if check_parameter correctly asserts when input parameters are missing."""
    out = ci.CheckInput(test, mode = mode)
    flag = out.check_parameter()
    assert flag == expected

@pytest.mark.parametrize(
    "test, mode, expected",
    [
        (dict(v=1), None, True),
        (dict(v=-1), None, False),
        (dict(v=1, mu=-1), None, False),
        (dict(v="no"), None, False),
        (dict(c_source=2), None, True),
        (dict(c_source="no"), None, False),
        (dict(c_source=np.array([[10, 2], [20, 3], [10, 2]])), None, True),
        (dict(c_source=np.array([[10, 2, 1], [20, 3, 1], [10, 2, 1]])), None, False),
        (dict(c_source=np.array([20, 10])), None, False),
        (dict(c_source=np.array([[-1, 2], [20, 3], [10, 2]])), None, False),
        (dict(R=1), None, True),
        (dict(R=0.5), None, False),
        (dict(n=0.5), None, True),
        (dict(n=-1), None, False),
        (dict(n=1.5), None, False),
    ])

def test_check_values(test, mode, expected):
    """Test if check_values correctly asserts when input parameters are of wrong type or value."""
    out = ci.CheckInput(test, mode = mode)
    flag = out.check_values()
    assert flag == expected
