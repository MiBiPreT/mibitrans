"""
Data class for handling field site specific transport parameters as input.

@author: Alraune Zech
"""

from dataclasses import dataclass
import logging

@dataclass
class DataClassTransport:
    # flow parameter (for advection)
    velo: float 
#    hydr_cond: float
#    hydr_grad: float    
    por: float = 1
    # dispersion parameters (for spreading)
    alpha_x: float = 1
    alpha_y: float = 1
    alpha_z: float = 1
    plume_length: float = 1 #???
    # adsorption and decay parameters
    ret: float = 1 # linear equilibrium adsorption coefficient --> retardation factor
    rho: float = 1 # density of water?  #[kg/m3]
    koc: float = 1 # in [m3/kg]
    foc: float = 1 #[-]
    decay_rate: float = 1 #[/d] decay rate
    half_life: float = 1 #[d] half life time
    # source dimensions
    source_depth: float = 1 #[m]
    source_conc: float = 1 #[g/m3]
    source_mass_total: float = 1 #[kg]   
    # electron acceptors
    delta_oxygen: float = 1 #[g/m3]
    delta_nitrate: float = 1 #[g/m3]
    ferrous_iron: float = 1 #[g/m3]
    delta_sulfate: float = 1 #[g/m3]
    methane: float = 1 #[g/m3]
    # spatial and temporal dimensions of the modell --> outsource to later
    l_model: float = 1 #[m]
    w_model: float = 1 #[m]
    t_model: float = 1 #[d]
    # units
    unit_time: str = 'd'
    unit_dist: str = 'm'
    unit_mass: str = 'kg'


def standardize_input_dict(dictionary: dict,
                     verbose: bool = True
                     ) -> dict:
    """Format and structure input dictionary into a standardized dictionary.

    Args:
        dictionary (dict): Input dictionary.

    Returns:
        dict: Dictionary following standardized format.
    """
    params = {}
    unknown_keys = []

    logging.info('==============================================================')
    logging.info('Running standardization of input parameter dictionary')
    logging.info('==============================================================')

    # Convert every input dictionary key by looping over all keys
    for key_input, value_input in dictionary.items():
        key_in_known_keys = False
        for key_params, key_known in standard_names.items():
            # Look if input key is listed as a possible name for a parameter
            if key_input in key_known:
                params[key_params] = value_input
                key_in_known_keys = True
                # If input key is recognized, no need to continue this loop and go to next input key.
                break
        if not key_in_known_keys:
            unknown_keys.append(key_input)

    if len(unknown_keys) > 0:

        logging.info("Keys not recognized and not included in parameter dictionary:\n", unknown_keys)
    else:
        logging.info("All keys were recognized.")

    return params

standard_names = {
    "velo" : ["velo", "v", "velocity", "vel", "flow velocity", "flow_velocity"],
    "hydr_cond" : ["hydr_cond","k", "conductivity","hydraulic conductivity",
                   "hydraulicconductivity","hydraulic-conductivity","hydraulic_conductivity"],
    "hydr_grad" : ["hydr_grad","i","dh", "deltah", "delta h","delta_h","delta-h",
                   "gradient", "hydraulic gradient", "hydraulic_gradient",
                   "hydraulicgradient", "hydraulic-gradient"],
    "por" : ["por", "n", "porosity", "por", "theta"],
    "alpha_x" : ["alpha_x", "a_x", "dispersivity_x", "alpha_l", "longitudinal_dispersivity"],
    "alpha_y" : ["alpha_y", "a_y", "dispersivity_y", "alpha_t", "transverse_dispersivity"],
    "alpha_z" : ["alpha_z", "a_z", "dispersivity_z", "alpha_v", "vertical_dispersivity"],
    "plume_length" : ["plume_length", "lp", "plumelength","plume-length", "plume length", "p_length", "pl"],
    "ret" : ["ret", "r", "ret_factor", "ret-factor", "ret factor", "retfactor",
             "retention_factor", "retention-factor", "retention factor", "retentionfactor"
             "ret_rate","ret-rate","ret rate","retrate"
             "retention_rate","retention-rate","retention rate","retentionrate",
             "adsorption_rate","adsorption-rate","adsorption rate","adsorptionrate",
             "adsorption_factor","adsorption-factor","adsorption factor","adsorptionfactor",
              "retardation_factor","retardation-factor","retardation factor","retardationfactor"
             ],
    "rho" : ["rho", "density", "bulk_density", "soil_density", "soil_bulk_density"],
    "koc" : ["koc", "partition_coefficient", "partition_coeff", "partition_coef"],
    "foc" : ["foc", "Foc", "organic_carbon", "fraction_organic_carbon", "fraction_organic"],
    "decay_rate" : ["decay_rate","mu", "Mu", "decay_coefficient",  "decay", "decay_coeff", "decay_coef"],
    "half_life" : ["half_life","t_half",  "t_1/2", "solute_half_life"],
    "l_model" : ["l_model", "model_l", "model_length", "length_model"],
    "w_model" : ["w_model", "model_w", "model_width", "width_model"],
    "t_model": ["t_model", "model_t", "time_model", "model_time", "simulation_t", "simulation_time",
                "t_end", "end_t", "end_time"],
    "source_depth" : ["source_depth","d_source", "source_thickness", "thickness_source"],
    "source_conc" : ["source_conc", "c_source", "concentration_source", "conc_source", "source_c", "source_concentration", "source_conc",
                  "source_data", "initial_conc", "c_initial", "conc_initial", "initial_concentration"],
    "source_mass_total" : ["source_mass_total","m_total", "total_m", "total_mass", "mass_total", "soluble_mass", "soluble_m",
                 "m_soluble", "source_mass"],
    "delta_oxygen" : ["dO", "dO2", "DO", "DO2", "delta_O", "delta_oxygen", "d_oxygen", "delta_O2", "do", "do2",
            "Do", "Do2", "O2", "oxygen", "Oxygen", "o2"],
    "delta_nitrate" : ["dNO3", "DNO3", "dN", "dNO", "delta_nitrate", "delta_NO3", "dno3", "Dno3", "d_nitrate", "nitrate",
              "Nitrate", "NO3", "no3"],
    "ferrous_iron" : ["Fe2", "Fe", "ferrous_iron", "fe", "fe2", "fe2+", "iron", "Iron", "Ferrous_Iron", "ferrous_iron"],
    "delta_sulfate" : ["dSO4", "dSO", "dS", "DSO4", "dSO42-", "sulfate", "delta_sulfate", "d_sulfate", "SO4", "so4", "SO42-",
              "delta_SO4"],
    "methane" : ["CH4", "ch4", "methane", "Methane", "CH", "ch"],
    "unit_time": ["unit_time"],
    "unit_dist": ["unit_dist"],
    "unit_mass": ["unit_mass"],
}
