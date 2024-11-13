"""Author: Jorrit Bakker.

File testing functionality of analytical_equation module.
"""

import numpy as np
import pytest
import anatrans.transport.analytical_equation as ana


@pytest.mark.parametrize(
    "test, mode, expect_raises",
    [
        (dict(v=1, lp=2, R=3, d_source=4, c_source=5, m_total=6), None, None),
        (dict(v=1, lp=2, R=3, c_source=5, m_total=6), None, ValueError),
        (dict(v=1, lp=2, R=3, d_source=4, m_total=6), None, ValueError),
        (dict(v=1, lp=2, R=3, d_source=4, c_source=5), None, ValueError),
        (dict(), None, False),
        (dict(k=0.5, i=1, n=1.5, lp=2, R=3, d_source=4, c_source=5, m_total=6), None, None),
        (dict(k=0.5, n=1.5, lp=2, R=3, d_source=4, c_source=5, m_total=6), None, ValueError),
        (dict(v=1, alpha_x=1.5, alpha_y=2, alpha_z=2.5, R=3, d_source=4, c_source=5, m_total=6), None, None),
        (dict(v=1, lp=2, rho=2.5, Koc=3, foc=3.5, n=3.75, d_source=4, c_source=5, m_total=6), None, None),
        (dict(k=0.5, i=1, n=1.5, lp=2, R=3, d_source=4, c_source=5, m_total=6, mu=7), "linear_decay", None),
        (dict(k=0.5, i=1, n=1.5, lp=2, R=3, d_source=4, c_source=5, m_total=6), "linear_decay", ValueError),
        (dict(k=0.5, i=1, n=1.5, lp=2, R=3, d_source=4, c_source=5, m_total=6, dO=7, dNO3=8, Fe2=9, dSO4=10, CH4=11),
         "instant_reaction", None),
        (dict(k=0.5, i=1, n=1.5, lp=2, R=3, d_source=4, c_source=5, m_total=6), "instant_reaction", ValueError),
        (dict(k=0.5, i=1, n=1.5, lp=2, R=3, d_source=5, m_total=6, dO=7, dNO3=8, Fe2=9, CH4=11),
         "instant_reaction", ValueError),
    ])

def test_transport_missingvalue(test, mode, expect_raises) -> None:
    """Test if check_parameter correctly asserts when input parameters are missing."""
    # if expect_raises is not None:
    #     with pytest.raises(expect_raises):
    #         result = patient_normalise(np.array(test))
    #         npt.assert_allclose(result, np.array(expected), rtol=1e-2, atol=1e-2)
    # else:
    #     result = patient_normalise(np.array(test))
    #     npt.assert_allclose(result, np.array(expected), rtol=1e-2, atol=1e-2)

    if expect_raises is not None:
        with pytest.raises(expect_raises):
            ana.Transport(test, mode=mode)
    else:
        ana.Transport(test, mode=mode)
#
# @pytest.mark.parametrize(
#     "test, mode, expected",
#     [
#         (dict(v=1), None, True),
#         (dict(v=-1), None, False),
#         (dict(v=1, mu=-1), None, False),
#         (dict(v="no"), None, False),
#         (dict(c_source=2), None, True),
#         (dict(c_source="no"), None, False),
#         (dict(c_source=np.array([[10, 2], [20, 3], [10, 2]])), None, True),
#         (dict(c_source=np.array([[10, 2, 1], [20, 3, 1], [10, 2, 1]])), None, False),
#         (dict(c_source=np.array([20, 10])), None, False),
#         (dict(c_source=np.array([[-1, 2], [20, 3], [10, 2]])), None, False),
#     ])
#
# def test_check_values(test, mode, expected):
#     """Test if check_values correctly asserts when input parameters are of wrong type or value."""
#     out = ci.CheckInput(test, mode = mode)
#     flag = out.check_values()
#     assert flag == expected
