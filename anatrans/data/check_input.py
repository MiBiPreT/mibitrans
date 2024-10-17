"""Author: Jorrit Bakker.

Module evaluating if a dictionary contains all required (correct) parameters for analysis
"""


class CheckInput:
    def __init__(self, dictionary, mode = None, verbose = True):
        self.dict = dictionary
        self.mode = mode
        self.verbose = verbose

    def check_parameter(self):
        """Check if all required parameters are present for a specific mode."""
        keys = self.dict.keys()
        self.missing_params = []

        if not("v" in keys or set(["k", "i", "n"]).issubset(keys)):
            self.missing_params.append("v or (k, i and n)")

        if not("lp" in keys or set(["alpha_x", "alpha_y", "alpha_z"]).issubset(keys)):
            self.missing_params.append("lp or (alpha_x and alpha_y and alpha_z)")

        if not("R" in keys or set(["rho", "Koc", "foc"]).issubset(keys)):
            self.missing_params.append("R or (rho, Koc and foc)")

        if not "d_source" in keys:
            self.missing_params.append("d_source")

        if not "m_total" in keys:
            self.missing_params.append("m_total")

        if self.mode == "linear_decay":
            if not("mu" in keys or "t_half" in keys):
                self.missing_params.append("mu or t_half")

        elif self.mode == "instant_reaction":
            if not(set(["dO", "dNO3", "Fe2", "dSO4", "CH4"]).issubset(keys)):
                self.missing_params.append("dO, dNO3, Fe2, dSO4, CH4")

        elif self.mode and self.mode != "no_decay":
            if self.verbose:
                print("Mode is not recognized. Only checking for standard requirements.")

        if self.verbose and self.missing_params:
            print("The following parameters are missing:", self.missing_params)

        if self.missing_params:
            success_flag = False
        else:
            success_flag = True

        return success_flag

    def check_values(self):
        print("b")