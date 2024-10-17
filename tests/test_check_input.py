"""Author: Jorrit Bakker.

Module handling data input in the form of a dictionary.
"""

import pytest
import anatrans.data.check_input as ci

@pytest.mark.parametrize(
    "test, mode, expected",
    [
        (dict(v = 1, lp = 2, R = 3, d_source = 4, m_total = 5), None, True),
        (dict(v = 1, lp = 2, R = 3, d_source = 4), None, False),
        (dict(k = 0.5, i = 1, n = 1.5, lp = 2, R = 3, d_source = 4, m_total = 5), None, True),
        (dict(k=0.5, n=1.5, lp=2, R=3, d_source=4, m_total=5), None, False),
        (dict(v = 1, alpha_x = 1.5, alpha_y = 2, alpha_z = 2.5, R=3, d_source=4, m_total=5), None, True),
        (dict(v=1, lp=2, rho = 2.5, Koc = 3, foc = 3.5, d_source=4, m_total=5), None, True),
        (dict(k=0.5, i=1, n=1.5, lp=2, R=3,d_source=4, m_total=5,  mu=6), "linear_decay", True),
        (dict(k=0.5, i=1, n=1.5, lp=2, R=3, d_source=4, m_total=5), "linear_decay", False),
        (dict(k=0.5, i=1, n=1.5, lp=2, R=3, d_source=5, m_total=6, dO = 7, dNO3 = 8, Fe2 = 9, dSO4 = 10, CH4 = 11)
         , "instant_reaction", True),
        (dict(k=0.5, i=1, n=1.5, lp=2, R=3, d_source=5, m_total=6), "instant_reaction", False),
    ])

def test_check_parameter(test, mode, expected) -> None:
    """Test if from_dict gives expected output for various input dictionaries."""
    out = ci.CheckInput(test, mode = mode)
    flag = out.check_parameter()
    assert flag == expected