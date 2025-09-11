"""Author: Jorrit Bakker.

Module calculating the solution to the Domenico (1987) analytical model adapted in BIOSCREEN, for different scenarios.
"""
import warnings
import numpy as np
from scipy.special import erf, erfc
from mibitrans.data.parameter_information import util_to_conc_name


class Domenico:
    """Parent class that for all analytical solutions using on the Domenico analytical model."""
    def __init__(self,
                 hydrological_parameters,
                 adsorption_parameters,
                 source_parameters,
                 model_parameters,
                 verbose = False
                 ):
        self.hyd_pars = hydrological_parameters
        self.ads_pars = adsorption_parameters
        self.src_pars = source_parameters
        self.mod_pars = model_parameters
        self.verbose = verbose

        # One-dimensional model domain arrays
        self.x = np.arange(0, self.mod_pars.model_length + self.mod_pars.dx, self.mod_pars.dx)
        self.y = self._calculate_y()
        self.t = np.arange(0, self.mod_pars.model_time + self.mod_pars.dt, self.mod_pars.dt)

        # Three-dimensional model domain arrays
        self.xxx = np.tile(self.x, (len(self.t), len(self.y), 1))
        self.yyy = np.tile(self.y[:, None], (len(self.t), 1, len(self.x)))
        self.ttt = np.tile(self.t[:, None, None], (1, len(self.y), len(self.x)))
        # cxyt is concentration output array
        self.cxyt = np.zeros(self.xxx.shape)

        # Calculate retardation (if not already specified in adsorption_parameters)
        self.ads_pars.calculate_retardation(self.hyd_pars.porosity)
        # Calculate retarded velocity
        self.rv = self.hyd_pars.velocity / self.ads_pars.retardation

        self.k_source = self.calculate_source_decay()
        self.y_source = self.src_pars.source_zone_boundary
        # Subtract outer source zones from inner source zones
        self.c_source = self.src_pars.source_zone_concentration.copy()
        self.c_source[:-1] = self.c_source[:-1] - self.c_source[1:]

    # Could be moved to general file with functions, since every solution will use this
    def _calculate_y(self):
        if self.mod_pars.model_width > 2 * self.src_pars.source_zone_boundary[-1]:
            y = np.arange(-self.mod_pars.model_width / 2, self.mod_pars.model_width / 2, self.mod_pars.dy)
        else:
            y = np.arange(-self.src_pars.source_zone_boundary[-1], self.src_pars.source_zone_boundary[-1] + self.mod_pars.dy, self.mod_pars.dy)
            warnings.warn("Source zone boundary is larger than model width. Model width adjusted to fit entire source zone.")
        return y

    def calculate_source_decay(self, biodegradation_capacity=0):
        if isinstance(self.src_pars.total_mass, (float, int)):
            y_src = np.zeros(len(self.src_pars.source_zone_boundary)+1)
            y_src[1:]= self.src_pars.source_zone_boundary
            c_src = self.src_pars.source_zone_concentration
            Q = self.hyd_pars.velocity * self.hyd_pars.porosity * self.src_pars.depth * np.max(y_src) * 2

            weighted_conc = np.zeros(len(self.src_pars.source_zone_boundary))
            for i in range(len(self.src_pars.source_zone_boundary)):
                weighted_conc[i] = (y_src[i+1] - y_src[i]) * c_src[i]

            c0_avg = biodegradation_capacity + np.sum(weighted_conc) / np.max(y_src)
            print("c0_avg", c0_avg)
            k_source = Q * c0_avg / self.src_pars.total_mass
        # If source mass is not a float, it is an infinite source, therefore, no source decay takes place.
        else:
            k_source = 0

        return k_source

    def _eq_x_term(self, decay_sqrt=1):
        return erfc((self.xxx - self.hyd_pars.velocity * self.ttt * decay_sqrt)
                    / (2 * np.sqrt(self.hyd_pars.alpha_x * self.hyd_pars.velocity * self.ttt)))

    def _eq_additional_x(self):
        return(np.exp(self.xxx * self.hyd_pars.velocity / (self.hyd_pars.alpha_x * self.hyd_pars.velocity))
                       * erfc(self.xxx + self.hyd_pars.velocity * self.ttt/ (2 * np.sqrt(self.hyd_pars.alpha_x * self.hyd_pars.velocity * self.ttt))))

    def _eq_z_term(self):
        inner_term = self.src_pars.depth / (2 * np.sqrt(self.hyd_pars.alpha_z * self.xxx))
        return erf(inner_term) - erf(-inner_term)

    def _eq_source_decay(self):
        term = np.exp(-self.k_source * (self.ttt - self.xxx / self.hyd_pars.velocity))
        # Term can be max 1; can not have 'generation' of solute ahead of advection.
        return np.where(term > 1, 1, term)

    def _eq_y_term(self, i):
        div_term = 2 * np.sqrt(self.hyd_pars.alpha_y * self.xxx)
        term = erf((self.yyy + self.y_source[i]) / div_term) - erf((self.yyy - self.y_source[i]) / div_term)
        term[np.isnan(term)] = 0
        return term


class no_decay(Domenico):
    def __init__(self,
                 hydrological_parameters,
                 adsorption_parameters,
                 source_parameters,
                 model_parameters,
                 verbose=False
                 ):
        super().__init__(hydrological_parameters, adsorption_parameters, source_parameters, model_parameters, verbose)

        self._calculate()

    def _calculate(self):
        with np.errstate(divide="ignore", invalid="ignore"):
            x_term = self._eq_x_term()
            additional_x = self._eq_additional_x()
            z_term = self._eq_z_term()
            source_decay = self._eq_source_decay()
            for i in range(len(self.c_source)):
                y_term = self._eq_y_term(i)
                cxyt = 1/8 * self.c_source[i] * source_decay * (x_term + additional_x) * y_term * z_term
                self.cxyt += cxyt


class linear_decay(Domenico):
    def __init__(self,
                 hydrological_parameters,
                 adsorption_parameters,
                 degradation_parameters,
                 source_parameters,
                 model_parameters,
                 verbose = False
                 ):
        super().__init__(hydrological_parameters, adsorption_parameters, source_parameters, model_parameters, verbose)
        self.deg_pars = degradation_parameters
        self.deg_pars.require_linear_decay()
        self._calculate()

    def _calculate(self):
        with np.errstate(divide="ignore", invalid="ignore"):
            decay_sqrt = np.sqrt(1 + 4 * self.deg_pars.decay_rate * self.hyd_pars.alpha_x / self.hyd_pars.velocity)
            decay_term = np.exp(self.xxx * (1 - decay_sqrt) / (self.hyd_pars.alpha_x * 2))
            x_term = self._eq_x_term(decay_sqrt)
            z_term = self._eq_z_term()
            source_decay = self._eq_source_decay()
            for i in range(len(self.c_source)):
                y_term = self._eq_y_term(i)
                cxyt = 1/8 * self.c_source[i] * source_decay * decay_term * x_term * y_term * z_term
                self.cxyt += cxyt


class instant_reaction(Domenico):
    def __init__(self,
                 hydrological_parameters,
                 adsorption_parameters,
                 degradation_parameters,
                 source_parameters,
                 model_parameters,
                 verbose = False
                 ):
        super().__init__(hydrological_parameters, adsorption_parameters, source_parameters, model_parameters, verbose)
        self.deg_pars = degradation_parameters
        self.deg_pars.require_electron_acceptor()
        self.biodegradation_capacity = self._calculate_biodegradation_capacity()
        self.k_source = self.calculate_source_decay(self.biodegradation_capacity)
        # Source decay calculation uses self.c_source, therefore, addition of biodegradation_capacity to outer source zone after calculation of k_source
        self.c_source[-1] += self.biodegradation_capacity
        self.cxyt_noBC = 0
        self._calculate()

    def _calculate_biodegradation_capacity(self):
        biodegradation_capacity = 0

        for key, item in self.deg_pars.electron_acceptor_utilization.items():
            biodegradation_capacity += self.deg_pars.__dict__[util_to_conc_name[key]] / item

        return biodegradation_capacity

    def _calculate(self):
        with np.errstate(divide="ignore", invalid="ignore"):
            x_term = self._eq_x_term()
            additional_x = self._eq_additional_x()
            z_term = self._eq_z_term()
            source_decay = self._eq_source_decay()

            for i in range(len(self.c_source)):
                y_term = self._eq_y_term(i)
                cxyt = 1 / 8 * self.c_source[i] * source_decay * (x_term + additional_x) * y_term * z_term
                self.cxyt += cxyt

            self.cxyt_noBC = self.cxyt.copy()
            self.cxyt -= self.biodegradation_capacity
            self.cxyt = np.where(self.cxyt < 0, 0, self.cxyt)
