"""Author: Jorrit Bakker.

Module calculating the mass balance based on base parameters.
"""

import numpy as np
import anatrans.transport.analytical_equation as eq
from anatrans.analysis.parameter_calculations import calculate_acceptor_utilization
from anatrans.analysis.parameter_calculations import calculate_biodegradation_capacity
from anatrans.analysis.parameter_calculations import calculate_source_decay
from anatrans.analysis.parameter_calculations import calculate_source_decay_instant


class MassBalance:
    """Calculate contaminant mass balance across model compartments."""
    def __init__(self,
                 parameters : dict,
                 mode: str = "no_decay",
                 dx : float = None,
                 dy : float = None,
                 dt : float = None,
                 ) -> None:
        """Initialise the class and internal variables.

        Args:
            parameters (dict) : Dictionary with transport parameters.
            mode (str) : Type of analytical model to be used. Default is no decay model.
            dx (float) : Model step size in x direction. If left empty,
            reasonable value will be calculated based on modeled area length. Default is None.
            dy (float) : Model step size in y direction. If left empty,
            reasonable value will be calculated based on modeled area width. Default is None.
            dt (float) : Model step size for time. If left empty,
            reasonable value will be calculated based on simulation time. Default is None.

        """
        self.pars = parameters
        self.mode = mode
        self.dx = dx
        self.dy = dy
        self.dt = dt
        self.mass_balance_dict = {
        }


    def balance(self, time = None) -> dict:
        """Calculates mass balance at a certain time point using the analytical equation for specified model type.

        Returns:
            mass_balance_dict (dict) : Dictionary containing mass for each mass balance component
            relevant to the model type.
        """
        obj_nodecay = eq.Transport(self.pars, dx = self.dx, dy = self.dy, dt = self.dt, mode = "no_decay")
        cxyt_nd, x, y, t = obj_nodecay.domenico()

        # If time point is specified, closest point in time array t is taken.
        # If not specified, defaults to last time point.
        if time is not None:
            time_pos = np.argmin(abs(t - time))
            self.mass_balance_dict["time"] = t[time_pos]
        else:
            time_pos = -1
            self.mass_balance_dict["time"] = t[time_pos]

        ksource = calculate_source_decay(self.pars)

        # Total source mass at t=0
        M_source_0 = self.pars["m_total"] * 1000
        self.mass_balance_dict["source_mass_0"] = M_source_0

        # Total source mass at t=t, for the no decay model
        M_source_t = M_source_0 * np.exp(-ksource * t[time_pos])
        self.mass_balance_dict["source_mass_t"] = M_source_t

        # Change in source mass at t=t, due to source decay by transport
        M_source_delta = M_source_0 - M_source_t
        self.mass_balance_dict["source_mass_change"] = M_source_delta

        # Volume of single cell, as dx * dy * source thickness
        cellsize = abs(x[1] - x[2]) * abs(y[1] - y[2]) * self.pars["d_source"]

        # Plume mass of no decay model; concentration is converted to mass by multiplying by cellsize and pore space.
        plume_mass_nodecay = np.sum(cxyt_nd[time_pos, :, 1:] * cellsize * self.pars["n"])
        self.mass_balance_dict["plume_mass_no_decay"] = plume_mass_nodecay

        # Difference between current plume mass and change in source mass must have been transported outside of model
        # extent for no decay scenarios; preservation of mass.
        if M_source_delta - plume_mass_nodecay < 0:
            transport_outside_extent_nodecay = 0
            self.mass_balance_dict["transport_outside_extent"] = transport_outside_extent_nodecay
        else:
            transport_outside_extent_nodecay = M_source_delta - plume_mass_nodecay
            self.mass_balance_dict["transport_outside_extent_nodecay"] = transport_outside_extent_nodecay

        if self.mode == "linear_decay":
            obj_decay = eq.Transport(self.pars, dx = self.dx, dy = self.dy, dt = self.dt, mode = "linear_decay")
            cxyt_dec, x, y, t = obj_decay.domenico()

            # Plume mass of linear decay model.
            plume_mass_lindecay = np.sum(cxyt_dec[time_pos, :, 1:] * cellsize * self.pars["n"])
            self.mass_balance_dict["plume_mass_linear_decay"] = plume_mass_lindecay

            # Calculate transport out of model extent linear decay as fraction of transport out of model for no decay
            # model, scaled by ratio between no decay and linear decay plume mass.
            transport_outside_extent_lindecay = (transport_outside_extent_nodecay * plume_mass_lindecay
                                                 / plume_mass_nodecay)
            self.mass_balance_dict["transport_outside_extent_lineardecay"] = transport_outside_extent_lindecay

            # Contaminant mass degraded by linear decay is diffrence plume mass no and linear decay plus difference in
            # mass transported outside model extent by no and linear decay.
            degraded_mass = (plume_mass_nodecay - plume_mass_lindecay + transport_outside_extent_nodecay
                               - transport_outside_extent_lindecay)
            self.mass_balance_dict["plume_mass_degraded_linear"] = degraded_mass

        elif self.mode == "instant_reaction":
            obj_inst = eq.Transport(self.pars, dx = self.dx, dy = self.dy, dt = self.dt, mode = "instant_reaction")
            cxyt_inst, x, y, t = obj_inst.domenico()

            # Matrix with concentration values before subtraction of biodegradation capacity
            cxyt_noBC = obj_inst.cxyt_noBC

            BC = calculate_biodegradation_capacity(self.pars)
            ksource_inst = calculate_source_decay_instant(self.pars, BC)

            # Total source mass at t=t, for the instant reaction model
            M_source_t_inst = M_source_0 * np.exp(-ksource_inst * t[time_pos])
            self.mass_balance_dict["source_mass_instant_t"] = M_source_t_inst

            # Change in source mass at t=t due to source decay by transport and by biodegradation
            M_source_delta = M_source_0 - M_source_t_inst
            self.mass_balance_dict["source_mass_instant_change"] = M_source_delta

            # Plume mass without biodegradation according to the instant degradation model
            plume_mass_nodecay = np.sum(cxyt_noBC[time_pos, :, 1:] * cellsize * self.pars["n"])
            self.mass_balance_dict["plume_mass_no_decay_instant_reaction"] = plume_mass_nodecay

            # Plume mass with biodegradation according to the instant degradation model
            plume_mass_inst = np.sum(cxyt_inst[time_pos, :, 1:] * cellsize * self.pars["n"])
            self.mass_balance_dict["plume_mass_instant_reaction"] = plume_mass_inst

            # Assumption: all mass difference between instant degradation model with biodegradation and
            # instant degradation model with biodegradation is caused by degradation.
            degraded_mass = plume_mass_nodecay - plume_mass_inst
            self.mass_balance_dict["plume_mass_degraded_instant"] = degraded_mass

            # Weight fraction of electron acceptor used for degradation and degraded contaminant
            mass_fraction_electron_acceptor = calculate_acceptor_utilization(self.pars)

            # Change in total mass of each electron acceptor
            electron_acceptor_mass_change = mass_fraction_electron_acceptor * degraded_mass
            self.mass_balance_dict["electron_acceptor_mass_change"] = electron_acceptor_mass_change

        return self.mass_balance_dict
