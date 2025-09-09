"""Author: Jorrit Bakker.

Module calculating the solution to the Domenico (1987) analytical model adapted in BIOSCREEN, for different scenarios.
"""

from mibitrans.data.read import HydrologicalParameters, AdsorptionDegradationParameters, ModelParameters

def no_decay(
    hydrological_parameters,
    adsorption_degradation_parameters,
    source_parameters,
    model_parameters = None,
):

    adsorption_degradation_parameters.calculate_retardation(hydrological_parameters.porosity)

    return None

def linear_decay(
    hydrological_parameters,
    adsorption_degradation_parameters,
    source_parameters,
    model_parameters,
):
    return None

def instant_reaction(
    hydrological_parameters,
    adsorption_degradation_parameters,
    source_parameters,
    model_parameters = None,
):
    return None

# def _evaluate_model_domain(
#     hydrological_parameters,
#     adsorption_degradation_parameters,
#     source_parameters,
#     model_parameters = None
# ):
