"""Author: Jorrit Bakker.

File containing valid parameter names for input in the form of a dictionary.
"""

fromdict_dictionary = {
    "v" : ["v", "V", "velocity", "vel", "velo"],
    "k" : ["k", "K", "hydraulic conductivity", "conductivity", "hydraulic_conductivity"],
    "i" : ["i", "I", "gradient", "hydraulic gradient", "hydraulic_gradient"],
    "n" : ["n", "N", "porosity", "por", "theta"],
    "alpha_x" : ["alpha_x", "a_x", "dispersivity_x", "alpha_l", "longitudinal_dispersivity"],
    "alpha_y" : ["alpha_y", "a_y", "dispersivity_y", "alpha_t", "transverse_dispersivity"],
    "alpha_z" : ["alpha_z", "a_z", "dispersivity_z", "alpha_v", "vertical_dispersivity"],
    "lp" : ["lp", "Lp", "plume_length", "plume length", "p_length", "pl", "Pl"],
    "R" : ["R", "r", "retention_factor", "adsorption_factor", "retention_rate", "adsorption_rate",
           "ret_factor", "ret_rate"],
    "rho" : ["rho", "Rho", "density", "bulk_density", "soil_density", "soil_bulk_density"],
    "Koc" : ["Koc", "koc", "partition_coefficient", "partition_coeff", "partition_coef"],
    "foc" : ["foc", "Foc", "organic_carbon", "fraction_organic_carbon", "fraction_organic"],
    "mu" : ["mu", "Mu", "decay_coefficient", "decay_rate", "decay", "decay_coeff", "decay_coef"],
    "t_half" : ["t_half", "half_life", "t_1/2", "solute_half_life"],
    "l_model" : ["l_model", "model_l", "model_length", "length_model"],
    "w_model" : ["w_model", "model_w", "model_width", "width_model"],
    "t_model": ["t_model", "model_t", "time_model", "model_time", "simulation_t", "simulation_time",
                "t_end", "end_t", "end_time"],
    "d_source" : ["d_source", "source_thickness", "thickness_source"],
    "m_total" : ["m_total", "total_m", "total_mass", "mass_total", "soluble_mass", "soluble_m", "m_soluble"],
    "dO" : ["dO", "dO2", "DO", "DO2", "delta_O", "delta_oxygen", "d_oxygen", "delta_O2", "do", "do2",
            "Do", "Do2", "O2", "oxygen", "Oxygen", "o2"],
    "dNO3" : ["dNO3", "DNO3", "dN", "dNO", "delta_nitrate", "delta_NO3", "dno3", "Dno3", "d_nitrate", "nitrate",
              "Nitrate", "NO3", "no3"],
    "Fe2+" : ["Fe2+", "Fe2", "Fe", "ferrous_iron", "fe", "fe2", "fe2+", "iron", "Iron", "Ferrous_Iron", "ferrous_iron"],
    "dSO4" : ["dSO4", "dSO", "dS", "DSO4", "dSO42-", "sulfate", "delta_sulfate", "d_sulfate", "SO4", "so4", "SO42-",
              "delta_SO4"],
    "CH4" : ["CH4", "ch4", "methane", "Methane", "CH", "ch"],
}

# Missing parameters: Modeled area length, modeled area width, simulation time, source thickness
# source zone concentrations, soluble mass

# Maybe only