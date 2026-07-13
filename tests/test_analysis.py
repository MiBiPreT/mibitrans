import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from mibitrans.analysis.differences import absolute_error
from mibitrans.analysis.differences import check_shape
from mibitrans.analysis.differences import comparison_plot
from mibitrans.analysis.differences import mask
from mibitrans.analysis.differences import mean_absolute_difference
from mibitrans.analysis.differences import mean_relative_difference
from mibitrans.analysis.differences import mean_sum_absolute_difference
from mibitrans.analysis.differences import relative_error
from mibitrans.analysis.differences import rmse
from mibitrans.analysis.differences import sensitivity_models
from mibitrans.transport.model_parent import Results


@pytest.mark.parametrize(
    "parameters, expected",
    [
        (
            dict(a=np.array([[2, 3], [4, 5]]), b=np.array([[3, 4], [5, 6]]), cutoff=4.5),
            (np.array([[np.nan, np.nan], [4, 5]]), np.array([[np.nan, np.nan], [5, 6]])),
        ),
        # If no cutoff given, should return input
        (
            dict(a=np.array([[0.01, 0.02], [-0.01, -0.02]]), b=np.array([[0.015, 0.01], [-0.03, -0.09]]), cutoff=None),
            (np.array([[0.01, 0.02], [-0.01, -0.02]]), np.array([[0.015, 0.01], [-0.03, -0.09]])),
        ),
    ],
)
def test_mask(parameters, expected) -> None:
    """Test if mask function correctly filters out values below cutoff value."""
    a, b = mask(**parameters)
    # By default, np.nan == np.nan gives False, specify that nan values are allowed
    assert a == pytest.approx(expected[0], nan_ok=True)
    assert b == pytest.approx(expected[1], nan_ok=True)


@pytest.mark.parametrize(
    "array_a, array_b, expected",
    [
        (np.array([[3, 2, 1], [3, 2, 1]]), np.array([[4, 3, 2], [4, 3, 2]]), None),
        (np.array([[3, 2, 1], [3, 2, 1]]), np.array([[4, 3, 2, 1], [4, 3, 2, 1]]), ValueError),
        (np.array([[3, 2, 1], [3, 2, 1]]), np.array([4, 3, 2, 1]), ValueError),
    ],
)
def test_shape(array_a, array_b, expected) -> None:
    """Test if array shapes are checked correctly."""
    if expected is not None:
        with pytest.raises(expected):
            check_shape(array_a, array_b)
    else:
        check_shape(array_a, array_b)


@pytest.mark.parametrize(
    "a, b, expected",
    [
        (np.array([[3, 2, 1], [3, 2, 1]]), np.array([[-4, 3, 5], [4, -3, 1]]), np.array([[7, 1, 4], [1, 5, 0]])),
        (4, -6, 10),
    ],
)
def test_absolute_error(a, b, expected) -> None:
    """Test if absolute error is evaluated correctly."""
    assert absolute_error(a, b) == pytest.approx(expected)


@pytest.mark.parametrize(
    "a, b, expected",
    [
        (
            np.array([[3, 2, 1], [3, 2, 1]]),
            np.array([[-4, 3, 5], [4, -3, 1]]),
            np.array([[2.33333333, 0.5, 4.0], [0.33333333, 2.5, 0.0]]),
        ),
        (4, 6, 0.5),
    ],
)
def test_relative_error(a, b, expected) -> None:
    """Test if absolute error is evaluated correctly."""
    assert relative_error(a, b) == pytest.approx(expected)


a_array = np.array([[[1, 2, 3], [4, 5, 6], [7, 8, 9]], [[10, 11, 12], [13, 14, 15], [16, 17, 18]]])
b_array = np.array(
    [[[1.5, 2.5, 3.5], [4.5, 5.6, 6.8], [7.1, 8.4, 9.3]], [[10.1, 11.2, 12.5], [13.3, 14.4, 15.9], [16.7, 17.7, 18]]]
)


@pytest.mark.parametrize(
    "a, b, axis, cutoff, expected",
    [
        (a_array, b_array, (0, 1, 2), None, 0.4444444444444444),
        (a_array, b_array, (0, 1), 3, np.array([0.34, 0.46, 0.5])),
    ],
)
def test_mean_absolute_difference(a, b, axis, cutoff, expected) -> None:
    """Test if mean absolute difference is calculated correctly."""
    assert mean_absolute_difference(a, b, axis, cutoff) == pytest.approx(expected)


@pytest.mark.parametrize(
    "a, b, axis, cutoff, expected",
    [
        (a_array, b_array, (0, 1, 2), None, 0.09216901970578441),
        (a_array, b_array, (0, 1), 3, np.array([0.04322253, 0.05158594, 0.0725])),
    ],
)
def test_mean_relative_difference(a, b, axis, cutoff, expected) -> None:
    """Test if mean absolute difference calculated correctly."""
    assert mean_relative_difference(a, b, axis, cutoff) == pytest.approx(expected)


@pytest.mark.parametrize(
    "a, b, cutoff, expected",
    [
        (a_array, b_array, 0, 1.3333333333333333),
        (a_array, b_array, 5, 1.4083333333333332),
    ],
)
def test_mean_sum_absolute_difference(a, b, cutoff, expected) -> None:
    """Test if function for sum of absolute differences works as expected."""
    assert mean_sum_absolute_difference(a, b, cutoff) == pytest.approx(expected)


@pytest.mark.parametrize(
    "a, b, axis, cutoff, expected",
    [
        (a_array, b_array, (0, 1, 2), 0, 0.5055250296034366),
        (a_array, b_array, (0, 1), 3, [0.41231056, 0.49193496, 0.58309519]),
    ],
)
def test_rmse(a, b, axis, cutoff, expected) -> None:
    """Test if RMSE is evaluated correctly."""
    assert rmse(a, b, axis, cutoff) == pytest.approx(expected)


@pytest.fixture(scope="module")
def test_sensitivty_models(
    test_hydro_pars,
    test_att_pars,
    test_source_pars,
    test_model_pars,
    x_dispersivity=[1, 0.5],
    y_dispersivity=[0.1, 0.05],
    z_dispersivity=[0.005, 0],
):
    """Output from sensitivity_models to check and use in testing plots."""
    list_mbt, list_ana, list_bio = sensitivity_models(
        test_hydro_pars,
        test_att_pars,
        test_source_pars,
        test_model_pars,
        x_dispersivity,
        y_dispersivity,
        z_dispersivity,
    )
    return list_mbt, list_ana, list_bio, x_dispersivity, y_dispersivity, z_dispersivity


def test_sensitivity_models_output(test_sensitivty_models):
    """Test if sensitivity_models function returns expected lists of model results."""
    list_mbt, list_ana, list_bio, _, _, _ = test_sensitivty_models
    assert isinstance(list_mbt, list)
    assert isinstance(list_mbt[0][0][0], Results)
    assert isinstance(list_ana, list)
    assert isinstance(list_ana[0][0][0], Results)
    assert isinstance(list_bio, list)
    assert isinstance(list_bio[0][0][0], Results)


def test_comparison_plot(test_sensitivty_models):
    """Test if comparison_plot creates axes objects and runs without error."""
    list_mbt, list_ana, _, x_dispersivity, y_dispersivity, z_dispersivity = test_sensitivty_models
    difference_methods = ["relative", "rmse", "absolute", "sum_mean"]
    for method in difference_methods:
        comparison_plot(
            list_mbt,
            list_ana,
            x_dispersivity,
            y_dispersivity,
            z_dispersivity,
            difference_method=method,
        )
        assert isinstance(plt.gca(), matplotlib.axes._axes.Axes)
        plt.clf()
