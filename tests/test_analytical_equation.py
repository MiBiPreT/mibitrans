"""Author: Jorrit Bakker.

File testing functionality of analytical_equation module.
"""

import numpy as np
import pytest
import anatrans.transport.analytical_equation as ana
from anatrans.data.parameter_information import testing_dictionary
from tests.test_example_results import testingdata_instantreaction
from tests.test_example_results import testingdata_lineardecay
from tests.test_example_results import testingdata_nodecay


@pytest.mark.parametrize(
    "input, stepsize",
    [
        (testing_dictionary, (10,5,1)),
        (testing_dictionary, (0, 0, 0)),
    ])

def test_transport_grid(input, stepsize) -> None:
    """Test if Transport class initial correctly makes the model grid."""
    dx, dy, dt = stepsize
    ini = ana.Transport(input, "no_decay", dx, dy, dt)

    source_y = input["c_source"][:, 0]

    if not isinstance(dx, (float, int)) or (dx <= 0):
        dx = input["l_model"] / 100
    if not isinstance(dx, (float, int)) or (dy <= 0):
        dy = input["w_model"] / 50
    if not isinstance(dx, (float, int)) or (dt <= 0):
        dt = input["t_model"] / 10

    x = np.arange(0, input["l_model"] + dx, dx)
    y = np.arange(-source_y[-1], source_y[-1] + dy, dy)
    t = np.arange(dt, input["t_model"] + dt, dt)

    test_xxx = np.tile(x, (len(t), len(y), 1))
    test_yyy = np.tile(y[:, None], (len(t), 1, len(x)))
    test_ttt = np.tile(t[:, None, None], (1, len(y), len(x)))
    test_cxyt = np.zeros(test_xxx.shape)
    print(test_xxx.shape, test_yyy.shape, test_ttt.shape)

    out_xxx = ini.xxx
    out_yyy = ini.yyy
    out_ttt = ini.ttt
    out_cxyt = ini.cxyt

    assert test_xxx == pytest.approx(out_xxx)
    assert test_yyy == pytest.approx(out_yyy)
    assert test_ttt == pytest.approx(out_ttt)
    assert test_cxyt == pytest.approx(out_cxyt)

@pytest.mark.parametrize(
    "input, stepsize, mode",
    [
        (testing_dictionary, (10,5,1), "no_decay"),
        (testing_dictionary, (10,5,1), "linear_decay"),
        (testing_dictionary, (10,5,1), "instant_reaction")
    ])

def test_transport_domenico(input, stepsize, mode) -> None:
    """Test if Transport class domenico model is correctly calculated."""
    dx, dy, dt = stepsize
    print(dx, dy, dt)
    obj = ana.Transport(testing_dictionary, mode, dx, dy, dt)
    out_cxyt, x, y, t = obj.domenico()
    if mode == "no_decay":
        test_cxyt = testingdata_nodecay
    elif mode == "linear_decay":
        test_cxyt = testingdata_lineardecay
    elif mode == "instant_reaction":
        print(out_cxyt)
        test_cxyt = testingdata_instantreaction
    else:
        test_cxyt = None

    assert test_cxyt == pytest.approx(out_cxyt)
