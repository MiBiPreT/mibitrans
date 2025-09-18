"""Author: Jorrit Bakker.

Module calculating the mass balance based on base parameters.
"""

import numpy as np
from mibitrans.analysis.parameter_calculations import calculate_utilization
from mibitrans.data.check_input import _check_model_type
from mibitrans.data.check_input import _time_check
from mibitrans.transport.domenico import Domenico
from mibitrans.transport.domenico import NoDecay


def mass_balance(model, time=None) -> dict:
    """Calculate contaminant mass balance across model compartments.

    Args:
        model (mibitrans.transport.domenico.Domenico) : Domenico transport model object.
        time (float) : Time at which to calculate mass balance. Default is the last time step.

    Returns:
        mass_balance_dict : Dictionary containing the mass balance elements of the given model.
    """
    _check_model_type(model, Domenico)
    time_pos = _time_check(model, time)
    mass_balance_dict = {}

    mass_balance_dict["time"] = model.t[time_pos]

    if hasattr(model, "deg_pars"):
        no_decay_model = NoDecay(model.hyd_pars, model.ads_pars, model.src_pars, model.mod_pars)
        if hasattr(model, "biodegradation_capacity"):
            mode = "instant_reaction"
        else:
            mode = "linear_decay"
    else:
        mode = "no_decay"
        no_decay_model = model

    # Total source mass at t=0
    M_source_0 = model.src_pars.total_mass
    mass_balance_dict["source_mass_0"] = M_source_0

    # Total source mass at t=t, for the no decay model
    M_source_t = M_source_0 * np.exp(-no_decay_model.k_source * model.t[time_pos])
    mass_balance_dict["source_mass_t"] = M_source_t

    # Change in source mass at t=t, due to source decay by transport
    M_source_delta = M_source_0 - M_source_t
    mass_balance_dict["source_mass_change"] = M_source_delta

    # Volume of single cell, as dx * dy * source thickness
    cellsize = abs(model.x[0] - model.x[1]) * abs(model.y[0] - model.y[1]) * model.src_pars.depth

    # Plume mass of no decay model; concentration is converted to mass by multiplying by cellsize and pore space.
    plume_mass_nodecay = np.sum(no_decay_model.cxyt[time_pos, :, 1:] * cellsize * model.hyd_pars.porosity)
    mass_balance_dict["plume_mass_no_decay"] = plume_mass_nodecay

    # Difference between current plume mass and change in source mass must have been transported outside of model
    # extent for no decay scenarios; preservation of mass.
    if M_source_delta - plume_mass_nodecay < 0:
        transport_outside_extent_nodecay = 0
        mass_balance_dict["transport_outside_extent"] = transport_outside_extent_nodecay
    else:
        transport_outside_extent_nodecay = M_source_delta - plume_mass_nodecay
        mass_balance_dict["transport_outside_extent_nodecay"] = transport_outside_extent_nodecay

    if mode == "linear_decay":
        # Plume mass of linear decay model.
        plume_mass_lindecay = np.sum(model.cxyt[time_pos, :, 1:] * cellsize * model.hyd_pars.porosity)
        mass_balance_dict["plume_mass_linear_decay"] = plume_mass_lindecay

        # Calculate transport out of model extent linear decay as fraction of transport out of model for no decay
        # model, scaled by ratio between no decay and linear decay plume mass.
        transport_outside_extent_lindecay = transport_outside_extent_nodecay * plume_mass_lindecay / plume_mass_nodecay
        mass_balance_dict["transport_outside_extent_lineardecay"] = transport_outside_extent_lindecay

        # Contaminant mass degraded by linear decay is diffrence plume mass no and linear decay plus difference in
        # mass transported outside model extent by no and linear decay.
        degraded_mass = (
            plume_mass_nodecay
            - plume_mass_lindecay
            + transport_outside_extent_nodecay
            - transport_outside_extent_lindecay
        )
        mass_balance_dict["plume_mass_degraded_linear"] = degraded_mass

    elif mode == "instant_reaction":
        # Total source mass at t=t, for the instant reaction model
        M_source_t_inst = M_source_0 * np.exp(-model.k_source * model.t[time_pos])
        mass_balance_dict["source_mass_instant_t"] = M_source_t_inst

        # Change in source mass at t=t due to source decay by transport and by biodegradation
        M_source_delta = M_source_0 - M_source_t_inst
        mass_balance_dict["source_mass_instant_change"] = M_source_delta

        # Plume mass without biodegradation according to the instant degradation model
        plume_mass_inst_nodecay = np.sum(model.cxyt_noBC[time_pos, :, 1:] * cellsize * model.hyd_pars.porosity)
        mass_balance_dict["plume_mass_no_decay_instant_reaction"] = plume_mass_inst_nodecay

        # Plume mass with biodegradation according to the instant degradation model
        plume_mass_inst = np.sum(model.cxyt[time_pos, :, 1:] * cellsize * model.hyd_pars.porosity)
        mass_balance_dict["plume_mass_instant_reaction"] = plume_mass_inst

        # Assumption: all mass difference between instant degradation model with biodegradation and
        # instant degradation model without biodegradation is caused by degradation.
        degraded_mass = plume_mass_inst_nodecay - plume_mass_inst
        mass_balance_dict["plume_mass_degraded_instant"] = degraded_mass

        # Weight fraction of electron acceptor used for degradation and degraded contaminant
        mass_fraction_electron_acceptor = calculate_utilization(model)

        # Change in total mass of each electron acceptor
        electron_acceptor_mass_change = mass_fraction_electron_acceptor * degraded_mass
        mass_balance_dict["electron_acceptor_mass_change"] = electron_acceptor_mass_change

    return mass_balance_dict
