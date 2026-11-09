"""Author: Jorrit Bakker.

Module containing various methods that takes a dictionary of parameters as input and calculates the proper values that
can be used in transport equations.
"""

import warnings
import numpy as np
from scipy.optimize import brentq
from scipy.special import erf


def calculate_biodegradation_capacity(electron_acceptors, utilization_factor):
    """Determine biodegradation capacity based on electron acceptor concentrations and utilization factor."""
    return np.sum(electron_acceptors.array / utilization_factor.array)


def calculate_discharge_and_average_source_zone_concentration(model):
    """Calculate discharge through source zone and average source zone concentration, returned in respective order."""
    if model.mode == "instant_reaction":
        bc = model._biodegradation_capacity
    else:
        bc = 0
    y_src = np.zeros(len(model.source_parameters.source_zone_boundary) + 1)
    y_src[1:] = model.source_parameters.source_zone_boundary
    c_src = model.source_parameters.source_zone_concentration
    Q = (
        model.hydrological_parameters.velocity
        * model.hydrological_parameters.porosity
        * model.source_parameters.depth
        * np.max(y_src)
        * 2
    )

    weighted_conc = np.zeros(len(model.source_parameters.source_zone_boundary))
    for i in range(len(model.source_parameters.source_zone_boundary)):
        weighted_conc[i] = (y_src[i + 1] - y_src[i]) * c_src[i]

    c0_avg = bc + np.sum(weighted_conc) / np.max(y_src)

    return Q, c0_avg


def _calculate_discharge_and_average_source_zone_concentration(
    hydrological_parameters, source_parameters, biodegradation_capacity=0
):
    """Calculate discharge through source zone and average source zone concentration, returned in respective order."""
    y_src = np.zeros(len(source_parameters.source_zone_boundary) + 1)
    y_src[1:] = source_parameters.source_zone_boundary
    c_src = source_parameters.source_zone_concentration
    Q = (
        hydrological_parameters.velocity
        * hydrological_parameters.porosity
        * source_parameters.depth
        * np.max(y_src)
        * 2
    )

    weighted_conc = np.zeros(len(source_parameters.source_zone_boundary))
    for i in range(len(source_parameters.source_zone_boundary)):
        weighted_conc[i] = (y_src[i + 1] - y_src[i]) * c_src[i]

    c0_avg = biodegradation_capacity + np.sum(weighted_conc) / np.max(y_src)

    return Q, c0_avg


def calculate_source_depletion(hydrological_parameters, source_parameters, biodegradation_capacity=0):
    """Calculate source depletion rate based on input parameters."""
    if source_parameters.total_mass != np.inf:
        Q, c0_avg = _calculate_discharge_and_average_source_zone_concentration(
            hydrological_parameters, source_parameters, biodegradation_capacity
        )
        k_source = Q * c0_avg / source_parameters.total_mass
    # If source mass is infinite, no source depletion is considered
    else:
        k_source = 0

    return k_source


def calculate_decay_ratios(decay_rates, mass_ratios):
    """Calculate decay ratios between contaminants in chain decay."""
    n = len(decay_rates)
    decay_ratios = np.full((n, n), np.nan)
    for i in range(n - 1):
        for j in range(i + 1, n):
            decay_ratios[i, j] = mass_ratios[i] * decay_rates[i] / (decay_rates[i] - decay_rates[j])
    return decay_ratios


def transform_chain_concentrations(concentration_iterable: list[np.ndarray], decay_ratios: np.ndarray, inverse=False):
    """Convert concentrations for chain decay into the transformed domain, or vice versa."""
    ni = len(concentration_iterable) - 1
    transformed_concentrations = [0] * (ni + 1)
    transformed_concentrations[0] = concentration_iterable[0]
    for k in range(ni):
        transform_total = 0
        for m in range(k + 1):
            loop_iterable = transformed_concentrations[m] if inverse else concentration_iterable[m]
            transform_total += loop_iterable * decay_ratios[m : k + 1, k + 1].prod()
        if inverse:
            transform_total *= -1
        transformed_concentrations[k + 1] = concentration_iterable[k + 1] + transform_total

    return transformed_concentrations


def calculate_plume_length(
    velocity, alpha_x, alpha_y, alpha_z, width, depth, decay_rate, c0, stoichiometric_ea_concentration
) -> float:
    """Estimate plume length based on flow, transport and degradation.

    Args:
        velocity (float): Flow velocity in the direction of the groundwater gradient.
        alpha_x (float): Longitudinal dispersivity.
        alpha_y (float): Horizontal transverse dispersivity.
        alpha_z (float): Vertical transverse dispersivity.
        width (float): Contaminant source width (from 0 to Y).
        depth (float): Contaminant source depth.
        decay_rate (float): First order (linear) decay coefficient.
        c0 (float): Contaminant source concentration (for multiple source zones, use average source zone concentration).
        stoichiometric_ea_concentration (float): Electron acceptor concentration stoichiometrically equal to
            electron donor. In case no fringe degradation is modelled, it is equal to 0.

    Returns:
        plume_length (float): Estimated plume length at model steady state.
    """
    max_length = 100 * width
    if stoichiometric_ea_concentration > 0:
        limit = stoichiometric_ea_concentration
    else:
        limit = c0 * 0.01

    # Equation used is t->infinity for Anatrans core-fringe degradation equation at plume center. Reminicent of equation
    # 23 of Gutierrez-Neri et al. (2009), but corrected for critiques of Hunkeler et al. (2010).
    def residual(L):
        """Calculate for which plume length equation is 0 at t->infinity."""
        degradation_term = np.exp((L / (2 * alpha_x)) * (1 - np.sqrt(1 + 4 * decay_rate * alpha_x / velocity))) + (
            stoichiometric_ea_concentration / c0
        )
        horizontal_transverse_term = erf(width / (2 * np.sqrt(alpha_y * L))) if alpha_y > 0 else 1
        vertical_transverse_term = erf(depth / (2 * np.sqrt(alpha_z * L))) if alpha_z > 0 else 1

        return c0 * degradation_term * horizontal_transverse_term * vertical_transverse_term - limit

    i = 0
    while residual(max_length) > 0:
        i += 1
        max_length *= 2
        if i > 100:
            warnings.warn("Failed to determine plume length. Returning nan.")
            return np.nan

    return brentq(residual, 0.001, max_length)
