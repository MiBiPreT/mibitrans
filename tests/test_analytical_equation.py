"""Author: Jorrit Bakker.

File testing functionality of analytical_equation module.
"""

import numpy as np
import pytest
import anatrans.transport.analytical_equation as ana

# base_dictionary = {
#     "v": 113.8 / 3.281,
#     "alpha_x": 13.3 / 3.281,
#     "alpha_y": 1.3 / 3.281,
#     "alpha_z": 0,
#     "R": 1,
#     "l_model": 320 / 3.281,
#     "w_model": 100 / 3.281,
#     "t_model": 6,
#     "d_source": 10 / 3.281,
#     "c_source": np.array([[0,13.68], [7 / 3.281,2.508], [37 / 3.281,0.057], [65 / 3.281,0]]),
#     "m_total": 2000,
#     "n" : 0.3,
#     "t_half" : 0.15
# }
# @pytest.mark.parametrize(
#     "test, mode, expect_raises",
#     [
#         (base_dictionary, None, None),
#         ({k: v for k, v in {**base_dictionary}.items() if k != 'v'}, None, ValueError),
#         (base_dictionary.copy().update({"v":-1}), None, ValueError),
#     ])
#
# def test_transport_missingvalue(test, mode, expect_raises) -> None:
#     """Test if check_parameter correctly asserts when input parameters are missing."""
#     # if expect_raises is not None:
#     #     with pytest.raises(expect_raises):
#     #         result = patient_normalise(np.array(test))
#     #         npt.assert_allclose(result, np.array(expected), rtol=1e-2, atol=1e-2)
#     # else:
#     #     result = patient_normalise(np.array(test))
#     #     npt.assert_allclose(result, np.array(expected), rtol=1e-2, atol=1e-2)
#
#     if expect_raises is not None:
#         with pytest.raises(expect_raises):
#             ana.Transport(test, mode, 1 / 3.281, 1 /3.281, 1)
#     else:
#         ana.Transport(test, mode, 1 / 3.281, 1 /3.281, 1)
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
