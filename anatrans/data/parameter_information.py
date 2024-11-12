"""Author: Jorrit Bakker.

File containing various dictionaries used for evaluation of names, value types and units of input data.
"""

# Dictionary with possible input variable names as value and the variable name used in the package as key.
key_dictionary = {
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
    "c_source" : ["c_source", "concentration_source", "conc_source", "source_c", "source_concentration", "source_conc",
                  "source_data", "initial_conc", "c_initial", "conc_initial", "initial_concentration"],
    "m_total" : ["m_total", "total_m", "total_mass", "mass_total", "soluble_mass", "soluble_m", "m_soluble"],
    "dO" : ["dO", "dO2", "DO", "DO2", "delta_O", "delta_oxygen", "d_oxygen", "delta_O2", "do", "do2",
            "Do", "Do2", "O2", "oxygen", "Oxygen", "o2"],
    "dNO3" : ["dNO3", "DNO3", "dN", "dNO", "delta_nitrate", "delta_NO3", "dno3", "Dno3", "d_nitrate", "nitrate",
              "Nitrate", "NO3", "no3"],
    "Fe2" : ["Fe2", "Fe", "ferrous_iron", "fe", "fe2", "fe2+", "iron", "Iron", "Ferrous_Iron", "ferrous_iron"],
    "dSO4" : ["dSO4", "dSO", "dS", "DSO4", "dSO42-", "sulfate", "delta_sulfate", "d_sulfate", "SO4", "so4", "SO42-",
              "delta_SO4"],
    "CH4" : ["CH4", "ch4", "methane", "Methane", "CH", "ch"],
}

# Dictionary representing which data types are allowed for each variable.
datatype_dictionary = {
    # Following parameters may be either floats or integers
    "int_float" : ["v", "k", "i", "n", "alpha_x", "alpha_y", "alpha_z", "lp", "R", "rho",
                      "Koc", "foc", "mu", "t_half", "l_model","w_model", "t_model", "d_source",
                      "m_total", "dO", "dNO3", "Fe2", "dSO4", "CH4"],
    # Following parameters may be floats, integers or arrays
    "float_array" : ["c_source"]
}

# Dictionary containing every input variable to serve as example.
example_dictionary = {
    "v": 1,
    "k": 1,
    "i": 1,
    "n": 1,
    "alpha_x": 1,
    "alpha_y": 1,
    "alpha_z": 1,
    "lp": 1,
    "R": 1,
    "rho": 1,
    "Koc": 1,
    "foc": 1,
    "mu": 1,
    "t_half": 1,
    "l_model": 1,
    "w_model": 1,
    "t_model": 1,
    "d_source": 1,
    "c_source": 1,
    "m_total": 1,
    "dO": 1,
    "dNO3": 1,
    "Fe2": 1,
    "dSO4": 1,
    "CH4": 1,
}

# unit_dictionary = {
#     "distance" : ["km", "m", "dm", "cm", "mm"],
#     "large_mass" : ["kg", "g"],
#     "small_mass" : ["g", "mg", "ug"],
#     "time" :
# }
