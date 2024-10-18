"""Author: Jorrit Bakker.

Module evaluating if a dictionary contains all required (correct) parameters for analysis
"""
import numpy as np

class CheckInput:
    """Evaluates if input dictionary contains all required parameters and if they have the correct data type."""
    def __init__(self, dictionary, mode = None, verbose = True) -> None:
        """Initialize class."""
        self.dict = dictionary
        self.mode = mode
        self.verbose = verbose
        self.wrong_type = []
        self.wrong_value = []
        self.missing_params = []

    def check_parameter(self) -> bool:
        """Check if all required parameters are present for a specific mode."""


        self.keys = self.dict.keys()

        if not("v" in self.keys or {"k", "i", "n"}.issubset(self.keys)):
            self.missing_params.append("v or (k, i and n)")

        if not("lp" in self.keys or {"alpha_x", "alpha_y", "alpha_z"}.issubset(self.keys)):
            self.missing_params.append("lp or (alpha_x and alpha_y and alpha_z)")

        if not("R" in self.keys or {"rho", "Koc", "foc"}.issubset(self.keys)):
            self.missing_params.append("R or (rho, Koc and foc)")

        if "d_source" not in self.keys:
            self.missing_params.append("d_source")

        if "c_source" not in self.keys:
            self.missing_params.append("c_source")

        if "m_total" not in self.keys:
            self.missing_params.append("m_total")

        if self.mode == "linear_decay":
            if not("mu" in self.keys or "t_half" in self.keys):
                self.missing_params.append("mu or t_half")

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

        for key, value in self.dict.items():
            if key in datatype_dictionary["int_float"]:
                if not isinstance(value, int) and not isinstance(value, float):
                    self.wrong_type.append(key)
                elif value < 0:
                    self.wrong_value.append(key)

            elif key in datatype_dictionary["float_array"]:
                if not isinstance(value, np.ndarray) and not isinstance(value, float) and not isinstance(value, int):
                    self.wrong_type.append(key)
                elif isinstance(value, int) or isinstance(value, float):
                    if value < 0:
                        self.wrong_value.append(key)
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