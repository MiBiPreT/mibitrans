"""Author: Jorrit Bakker.

Module containing various methods that takes a dictionary of parameters as input and calculates the proper values that
can be used in transport equations.
"""
import numpy as np


def calculate_retardation(pars):
    """Give retardation factor depending on input parameters."""
    if "R" in pars.keys():
        r = pars["R"]
    else:
        r = 1 + (pars["rho"] / pars["n"]) * pars["Koc"] * pars["foc"]
    return r

def calculate_flow_velocity(pars):
    """Give flow velocity depending on input parameters."""
    if "v" in pars.keys():
        v = pars["v"]
    else:
        v = (pars["k"] * pars["i"]) / (pars["n"])
    return v

def calculate_dispersivity(pars):
    """Give dispersivity in each direction depending on input parameters."""
    if "alpha_x" in pars.keys():
        ax = pars["alpha_x"]
        ay = pars["alpha_y"]
        az = pars["alpha_z"]
    else:
        # Implement conversion from plume length to alpha x, y and z
        ax = 0.83 * np.log10(pars["lp"])**2.414
        ay = ax / 10
        az = ax / 100

    if ax < 1e-10:
        ax = 1e-10
    if ay < 1e-10:
        ay = 1e-10
    if az < 1e-10:
        az = 1e-10

    return ax, ay, az

def calculate_linear_decay(pars):
    """Give 1st order decay coefficient depending on input parameters."""
    if "mu" in pars.keys():
        mu = pars["mu"]
    else:
        mu = np.log(2) / pars["t_half"]
    return mu

def calculate_source_decay(pars):
    """Function that calculates the source zone decay constant."""
    source_y = pars["c_source"][:, 0]
    source_c = pars["c_source"][:, 1]
    Q = pars["v"] * pars["n"] * pars["d_source"] * np.max(source_y) * 2 * 1000

    C0_avg = 0
    for i in range(len(source_y) - 1):
        if i == 0:
            yc = source_y[i + 1] * source_c[i] * 2
        else:
            yc = (source_y[i+1] - source_y[i]) * source_c[i] * 2
        C0_avg += yc
    C0_avg = C0_avg / (np.max(source_y) * 2)

    k_source = Q * C0_avg / (pars["m_total"] * 1e6)

    return(k_source)
