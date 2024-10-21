"""Author: Jorrit Bakker.

Module evaluating if a dictionary contains all required (correct) parameters for analysis
"""
import numpy as np


class CheckInput:
    """Evaluates if input dictionary contains all required parameters and if they have the correct data type."""
    def __init__(self, dictionary, mode = None, verbose = True) -> None:
        """Initialize parameters.

        Args:
            dictionary (dict): Input dictionary
            mode (str, optional): Method of analysis that will be performed on dictionary. Defaults to None.
            verbose (bool, optional): Verbose mode. Defaults to True.
        """
        self.dict = dictionary
        self.keys = self.dict.keys()
        self.mode = mode
        self.verbose = verbose
        self.wrong_type = []
        self.wrong_value = []
        self.missing_params = []

    def check_parameter(self) -> bool:
        """Check if all required parameters are present for a specific mode."""
        # Either the flow velocity or hydraulic conductivity, hydraulic gradient and porosity are required by the model.
        if not("v" in self.keys or {"k", "i", "n"}.issubset(self.keys)):
            self.missing_params.append("v or (k, i and n)")

        # Either plume length or dispersivity in all directions are required by the model.
        if not("lp" in self.keys or {"alpha_x", "alpha_y", "alpha_z"}.issubset(self.keys)):
            self.missing_params.append("lp or (alpha_x and alpha_y and alpha_z)")

        # Either retardation factor or bulk density, partition coefficient and fraction organic carbon
        # are required by the model.
        if not("R" in self.keys or {"rho", "Koc", "foc"}.issubset(self.keys)):
            self.missing_params.append("R or (rho, Koc and foc)")

        # Source thickness is required by the model.
        if "d_source" not in self.keys:
            self.missing_params.append("d_source")

        # Initial concentration is required by the model.
        if "c_source" not in self.keys:
            self.missing_params.append("c_source")

        # Total soluble mass is required by the model.
        if "m_total" not in self.keys:
            self.missing_params.append("m_total")

        # Either decay coefficient or solute half-life is required for the linear decay model.
        if self.mode == "linear_decay":
            if not("mu" in self.keys or "t_half" in self.keys):
                self.missing_params.append("mu or t_half")

        # Electron acceptor and donor concentrations are required for the instant reaction model.
        elif self.mode == "instant_reaction":
            if not({"dO", "dNO3", "Fe2", "dSO4", "CH4"}.issubset(self.keys)):
                self.missing_params.append("dO, dNO3, Fe2, dSO4, CH4")

        elif self.mode and self.mode != "no_decay":
            if self.verbose:
                print("Mode is not recognized. Only checking for standard requirements.")

        if self.missing_params:
            success_flag = False
            if self.verbose:
                print("The following parameters are missing:", self.missing_params)
        else:
            success_flag = True
            if self.verbose:
                print("All required parameters are present.")

        return success_flag

    def check_values(self):
        """Check if value and value types are as expected."""
        from anatrans.data.parameter_information import datatype_dictionary

        # To prevent incorrect datatypes, input data types are compared to allowed data types for each parameter.
        for key, value in self.dict.items():
            if key in datatype_dictionary["int_float"]:
                if not isinstance(value, int) and not isinstance(value, float):
                    self.wrong_type.append(key)
                # Parameters designated as float or int can not have negative values
                elif value < 0:
                    self.wrong_value.append(key)

            elif key in datatype_dictionary["float_array"]:
                if not isinstance(value, np.ndarray) and not isinstance(value, float) and not isinstance(value, int):
                    self.wrong_type.append(key)

                elif isinstance(value, int) or isinstance(value, float):
                    if value < 0:
                        self.wrong_value.append(key)

                # If input is an array, it needs to be 2d and have 2 columns.
                elif isinstance(value, np.ndarray):
                    if np.min(value) < 0:
                        self.wrong_value.append(key)
                    elif len(value.shape) != 2:
                        self.wrong_value.append(key)
                    elif value.shape[1] != 2:
                        self.wrong_value.append(key)

        if self.wrong_type or self.wrong_value:
            success_flag = False
            if self.verbose:
                print("The following parameters are of the wrong type:", self.wrong_type)
                print("The following parameters have a value that is smaller than 0:", self.wrong_value)

        else:
            success_flag = True
            if self.verbose:
                print("All parameters are of the correct type and value.")
        return(success_flag)

    def check_units(self):
        print("c")