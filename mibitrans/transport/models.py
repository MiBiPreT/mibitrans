import warnings
import numpy as np
from scipy.integrate import quad
from scipy.integrate import quad_vec
from scipy.special import erf
from scipy.special import erfc
from scipy.special import erfcx
from mibitrans.data.check_input import validate_input_values
from mibitrans.transport.model_parent import Results
from mibitrans.transport.model_parent import Transport3D


class Mibitrans(Transport3D):
    """Model class using an exact analytical solution as described in Karanovic (2007), based on Wexler (1992).

    Karanovic (2007) implemented the Wexler (1992) exact analytical solution in the Excel based BIOSCREEN-AT, and added
    source depletion, akin to that implemented in its predecessor BIOSCREEN by Newell et al. (1997). The Mibitrans model
    allows for the same method as used in BIOSCREEN-AT, but expands it by allowing multiple source zones (by means of
    superposition) and including the instant reaction model. These were present in the original BIOSCREEN, but not
    reimplemented in BIOSCREEN-AT. Using a single source zone in this model, and not using the instant reaction option
    will make the Mibitrans solution resolve to the equation described in Karanovic (2007). Which in turn resolves to
    the Wexler (1992) solution if source depletion is disabled.

    Karanovic, M., Neville, C. J., & Andrews, C. B. (2007). BIOSCREEN‐AT: BIOSCREEN with an exact analytical solution.
    Groundwater, 45(2), 242-245.

    Newell, C. J., McLeod, R. K., & Gonzales, J. R. (1997). BIOSCREEN natural attenuation decision support
    system version 1.4 revisions, Tech. rep., U.S. EPA.

    Wexler, E. J. (1992). Analytical solutions for one-, two-, and three-dimensional solute transport in ground-water
    systems with uniform flow. US Government Printing Office.
    """

    def __init__(
        self,
        hydrological_parameters,
        attenuation_parameters,
        source_parameters,
        model_parameters,
        verbose=False,
    ):
        """Initialize model object.

        Args:
            hydrological_parameters (mibitrans.data.parameters.HydrologicalParameters) : Dataclass object containing
                hydrological parameters from HydrologicalParameters.
            attenuation_parameters (mibitrans.data.read.AttenuationParameters) : Dataclass object containing adsorption,
                degradation and diffusion parameters from AttenuationParameters.
            source_parameters (mibitrans.data.read.SourceParameters) : Dataclass object containing source parameters
                from SourceParameters.
            model_parameters (mibitrans.data.read.ModelParameters) : Dataclass object containing model parameters from
                ModelParameters.
            verbose (bool, optional): Verbose mode. Defaults to False.

        Attributes:
            mode (str) : Current model mode. Is 'linear' by default. Once instant reaction parameters are provided. Use
                this attribute to switch between 'linear' and 'instant_reaction' models.
            cxyt (np.ndarray) : Output array containing concentrations in model domain, in [g/m^3]. Indexed as [t,y,x]
            relative_cxyt (np.ndarray) : Output array with concentrations in model domain, divided by the maximum source
                zone concentration at t=0. Indexed as [t,y,x].
            x (np.ndarray) : Discretized model x-dimension, in [m].
            y (np.ndarray) : Discretized model y-dimension, in [y].
            t (np.ndarray) : Discretized model t-dimension, in [days].
            c_source (np.ndarray) : Nett source zone concentrations, accounting for source superposition, in [g/m^3].
            vr (float) : Retarded groundwater flow velocity, in [m/d].
            k_source (float) : Source zone decay rate, in [1/days].
            biodegradation_capacity (float) : Maximum capacity of biodegradation given provided electron acceptor
            concentrations, in [g/m^3].

        Methods:
            run : Run model with current parameters, returns Results object.
            sample : Calculate concentration at any given position and point in time.
            instant_reaction : Activate the instant reaction model by providing electron acceptor concentrations. And
                optionally electron acceptor utilization factors. Switch between model modes by using the mode
                attribute.

        Raises:
            TypeError : If input is not of the correct Dataclass.

        """
        super().__init__(hydrological_parameters, attenuation_parameters, source_parameters, model_parameters, verbose)

    @property
    def short_description(self):
        """Return short description of model type."""
        if self.biodegradation_capacity:
            return "Mibitrans Instant Reaction"
        else:
            return "Mibitrans Linear"

    def run(self):
        """Calculate the concentration for all discretized x, y and t using the analytical transport model."""
        self._check_model_mode_before_run()
        with np.errstate(divide="ignore", invalid="ignore"):
            self.cxyt = self._calculate_concentration_for_all_xyt()
        return Results(self)

    def sample(self, x_position, y_position, time):
        """Give concentration at any given position and point in time.

        Args:
            x_position (float): x position in domain extent [m].
            y_position (float): y position in domain extent [m].
            time (float): time for which concentration is sampled [days].

        Returns:
            concentration (float): concentration at given position and point in time [g/m^3].

        """
        # Different sample method than parent class, as field-wide calculations use array indices
        for par, value in locals().items():
            if par != "self":
                validate_input_values(par, value)

        self._pre_run_initialization_parameters()

        def integrand(t, sz):
            div_term = 2 * np.sqrt(self.disp_y * t**4)
            inner_term = self._src_pars.depth / (2 * np.sqrt(self.disp_z * t**4))
            integrand_results = (
                1
                / (t**3)
                * (
                    np.exp(
                        (-self.k_source - self._decay_rate) * t**4
                        - (x_position - self.rv * t**4) ** 2 / (4 * self.disp_x * t**4)
                    )
                    * (
                        erfc((y_position - self.y_source[sz]) / div_term)
                        - erfc((y_position + self.y_source[sz]) / div_term)
                    )
                    * (erfc(-inner_term) - erfc(inner_term))
                )
            )
            return integrand_results

        conc_array = np.zeros(len(self.c_source))
        error_array = np.zeros(len(self.c_source))
        time = time ** (1 / 4)
        with np.errstate(divide="ignore", invalid="ignore"):
            for sz in range(len(self.c_source)):
                integral_term, error = quad(integrand, 0, time, limit=10000, args=(sz,))
                source_term = (
                    self.c_source[sz] * x_position / (8 * np.sqrt(np.pi * self.disp_x)) * np.exp(-self.k_source * time)
                )
                conc_array[sz] = 4 * integral_term * source_term
                error_array[sz] = error
            concentration = np.sum(conc_array)
            if self._mode == "instant_reaction":
                concentration -= self.biodegradation_capacity
                if concentration < 0:
                    concentration = 0
        return concentration

    def _pre_run_initialization_parameters(self):
        super()._pre_run_initialization_parameters()
        self.disp_x = self._hyd_pars.alpha_x * self.rv + self._att_pars.diffusion
        self.disp_y = self._hyd_pars.alpha_y * self.rv + self._att_pars.diffusion
        self.disp_z = self._hyd_pars.alpha_z * self.rv + self._att_pars.diffusion
        # self.integral_term = np.zeros(self.ttt.shape)
        # Stores integral error for each time step and source zone
        self.error_size = np.zeros((len(self._src_pars.source_zone_boundary), len(self.t)))

    def _calculate_concentration_for_all_xyt(self):
        cxyt = self.cxyt.copy()
        for sz in range(len(self.c_source)):
            integral_sum = self._equation_term_integral(sz)
            source_term = self._equation_term_source(sz)
            cxyt[:, :, 1:] += integral_sum[:, :, 1:] * source_term
            # If x=0, equation resolves to c=0, therefore, x=0 needs to be evaluated separately
            cxyt[:, :, 0] += self._equation_term_source_x_is_zero(sz)[:, :, 0]
        if self._mode == "instant_reaction":
            self.cxyt_noBC = cxyt.copy()
            cxyt -= self.biodegradation_capacity
            cxyt = np.where(cxyt < 0, 0, cxyt)
        return cxyt

    def _equation_term_integral(self, sz):
        integral_term = np.zeros(self.cxyt.shape)
        for j in range(len(self.t)):
            if self.verbose:
                print("integrating for source zone ", sz, " and t =", self.t[j], "days")
            if j == 0:
                lower_bound = 0
            else:
                lower_bound = self.t[j - 1]
            upper_bound = self.t[j]
            integral_term[j, :, 1:], self.error_size[sz, j] = quad_vec(
                self._equation_integrand, lower_bound, upper_bound, limit=10000 // len(self.t), args=(sz,)
            )
        integral_sum = np.cumsum(integral_term, axis=0)
        return integral_sum

    def _equation_integrand(self, t, sz):
        term = 1 / (t ** (3 / 2)) * self._equation_term_x(t) * self._equation_term_y(t, sz) * self._equation_term_z(t)
        term[np.isnan(term)] = 0
        return term

    def _equation_term_x(self, t):
        term = np.exp(
            (-self.k_source - self._decay_rate) * t - (self.xxx[:, :, 1:] - self.rv * t) ** 2 / (4 * self.disp_x * t)
        )
        term[np.isnan(term)] = 0
        return term

    def _equation_term_y(self, t, sz):
        div_term = 2 * np.sqrt(self.disp_y * t)
        term = erfc((self.yyy - self.y_source[sz]) / div_term) - erfc((self.yyy + self.y_source[sz]) / div_term)
        term[np.isnan(term)] = 0
        return term

    def _equation_term_z(self, t):
        if t == 0 or self.disp_z == 0:
            inner_term = 2
        else:
            inner_term = self._src_pars.depth / (2 * np.sqrt(self.disp_z * t))
        return erfc(-inner_term) - erfc(inner_term)

    def _equation_term_source(self, sz):
        return (
            self.c_source[sz]
            * self.xxx[:, :, 1:]
            / (8 * np.sqrt(np.pi * self.disp_x))
            * np.exp(-self.k_source * self.ttt)
        )

    def _equation_term_source_x_is_zero(self, sz):
        # Select y-positions of current source zone
        zone_location = np.where(abs(self.yyy) <= self.y_source[sz], 1, 0)
        return self.c_source[sz] * zone_location * np.exp(-self.k_source * self.ttt)


class Anatrans(Transport3D):
    """Model class using an analytical solution based on Bear (1979), Domenico (1987) & Newell et al. (1997).

    Under the assumption that C(x,y,z,t) = C(x,t) * C(y,t) * C(z,t), the 3D ADE can be broken up in three separate
    differential equations which can be solved individually. For C(x,t) the solution is given in Bear (1979), C(y,t) and
    C(z,t) can be derived from Crank (1975). The equation used for Anatrans is the combination of these solutions, with
    addition of source depletion, source superposition and instant reaction model, described in Newell et al. (1997) and
    implemented in the BIOSCREEN screening model. The solution of Newell et al. (1997) is based of the Domenico (1987)
    solution, a truncated version of the equation described above, which introduces an error with a size dependent on
    the ratio of flow velocity and longitudinal dispersivity. Anatrans instead uses the fully untruncated version.

    Bear, J. 1979. Hydraulics of Ground Water. New York: McGraw-Hill.

    Crank, J. 1975. The mathematics of Diffusion. New York: Oxford University Press.

    Domenico, P. A. (1987). An analytical model for multidimensional transport of a decaying contaminant species.
    Journal of Hydrology, 91(1-2), 49-58.

    Newell, C. J., McLeod, R. K., & Gonzales, J. R. (1997). BIOSCREEN natural attenuation decision support
    system version 1.4 revisions, Tech. rep., U.S. EPA.
    """

    def __init__(
        self,
        hydrological_parameters,
        attenuation_parameters,
        source_parameters,
        model_parameters,
        verbose=False,
    ):
        """Initialize model object.

        Args:
            hydrological_parameters (mibitrans.data.parameters.HydrologicalParameters) : Dataclass object containing
                hydrological parameters from HydrologicalParameters.
            attenuation_parameters (mibitrans.data.read.AttenuationParameters) : Dataclass object containing adsorption,
                degradation and diffusion parameters from AttenuationParameters.
            source_parameters (mibitrans.data.read.SourceParameters) : Dataclass object containing source parameters
                from SourceParameters.
            model_parameters (mibitrans.data.read.ModelParameters) : Dataclass object containing model parameters from
                ModelParameters.
            verbose (bool, optional): Verbose mode. Defaults to False.

        Attributes:
            mode (str) : Current model mode. Is 'linear' by default. Once instant reaction parameters are provided. Use
                this attribute to switch between 'linear' and 'instant_reaction' models.
            cxyt (np.ndarray) : Output array containing concentrations in model domain, in [g/m^3]. Indexed as [t,y,x]
            relative_cxyt (np.ndarray) : Output array with concentrations in model domain, divided by the maximum source
                zone concentration at t=0. Indexed as [t,y,x].
            x (np.ndarray) : Discretized model x-dimension, in [m].
            y (np.ndarray) : Discretized model y-dimension, in [y].
            t (np.ndarray) : Discretized model t-dimension, in [days].
            c_source (np.ndarray) : Nett source zone concentrations, accounting for source superposition, in [g/m^3].
            vr (float) : Retarded groundwater flow velocity, in [m/d].
            k_source (float) : Source zone decay rate, in [1/days].
            biodegradation_capacity (float) : Maximum capacity of biodegradation given provided electron acceptor
            concentrations, in [g/m^3].

        Methods:
            run : Run model with current parameters, returns Results object.
            sample : Calculate concentration at any given position and point in time.
            instant_reaction : Activate the instant reaction model by providing electron acceptor concentrations. And
                optionally electron acceptor utilization factors. Switch between model modes by using the mode
                attribute.
            centerline : Plot center of contaminant plume of this model, at a specified time and y position.
            transverse : Plot concentration distribution as a line horizontal transverse to the plume extent.
            breakthrough : Plot contaminant breakthrough curve at given x and y position in model domain.
            plume_2d : Plot contaminant plume as a 2D colormesh, at a specified time.
            plume_3d : Plot contaminant plume as a 3D surface, at a specified time.

        Raises:
            TypeError : If input is not of the correct Dataclass.

        """
        super().__init__(hydrological_parameters, attenuation_parameters, source_parameters, model_parameters, verbose)
        if self._att_pars.diffusion != 0:
            warnings.warn("Domenico model does not consider molecular diffusion.", UserWarning)

    @property
    def short_description(self):
        """Short description of model type."""
        return "Anatrans model"

    def run(self):
        """Calculate the concentration for all discretized x, y and t using the analytical transport model."""
        self._check_model_mode_before_run()
        with np.errstate(divide="ignore", invalid="ignore"):
            self.cxyt = self._calculate_concentration_for_all_xyt(self.xxx, self.yyy, self.ttt)
        return Results(self)

    def sample(self, x_position, y_position, time):
        """Give concentration at any given position and point in time.

        Args:
            x_position (float): x position in domain extent [m].
            y_position (float): y position in domain extent [m].
            time (float): time for which concentration is sampled [days].

        Returns:
            concentration (float): concentration at given position and point in time [g/m^3].

        """
        for par, value in locals().items():
            if par != "self":
                validate_input_values(par, value)

        self._pre_run_initialization_parameters()

        if self.mode == "instant_reaction":
            save_c_noBC = self.cxyt_noBC.copy()
        x = np.array([x_position])
        y = np.array([y_position])
        t = np.array([time])
        concentration = self._calculate_concentration_for_all_xyt(x, y, t)[0]
        if self.mode == "instant_reaction":
            self.cxyt_noBC = save_c_noBC
        return concentration

    def _equation_term_x(self, xxx, ttt, decay_sqrt):
        return np.exp(xxx * (1 - decay_sqrt) / (self._hyd_pars.alpha_x * 2)) * erfc(
            (xxx - self.rv * ttt * decay_sqrt) / (2 * np.sqrt(self._hyd_pars.alpha_x * self.rv * ttt))
        )

    def _equation_term_additional_x(self, xxx, ttt, decay_sqrt):
        erfc_inner = (xxx + decay_sqrt * self.rv * ttt) / (2 * np.sqrt(self._hyd_pars.alpha_x * self.rv * ttt))
        # Additional term is prone to overflow of exp and underflow of erfc under certain parameter combinations.
        # To decrease cases, used erfcx. Where erfcx(a) = exp(a**2)*erfc(a) -> exp(b)*erfc(a) = exp(b - a**2) * erfcx(a)
        term = np.exp(xxx * (1 + decay_sqrt) * (1 / 2) / self._hyd_pars.alpha_x - erfc_inner**2) * erfcx(erfc_inner)
        return term

    def _equation_term_z(self, xxx):
        inner_term = self._src_pars.depth / (2 * np.sqrt(self._hyd_pars.alpha_z * xxx))
        return erf(inner_term) - erf(-inner_term)

    def _equation_term_source_decay(self, xxx, ttt):
        term = np.exp(-self.k_source * (ttt - xxx / self.rv))
        # Term can be max 1; can not have 'generation' of solute ahead of advection.
        return np.where(term > 1, 1, term)

    def _equation_term_y(self, i, xxx, yyy):
        div_term = 2 * np.sqrt(self._hyd_pars.alpha_y * xxx)
        term = erf((yyy + self.y_source[i]) / div_term) - erf((yyy - self.y_source[i]) / div_term)
        term[np.isnan(term)] = 0
        return term

    def _calculate_concentration_for_all_xyt(self, xxx, yyy, ttt):
        cxyt = 0
        decay_sqrt = np.sqrt(1 + 4 * self._decay_rate * self._hyd_pars.alpha_x / self.rv)
        x_term = self._equation_term_x(xxx, ttt, decay_sqrt)
        additional_x = self._equation_term_additional_x(xxx, ttt, decay_sqrt)
        z_term = self._equation_term_z(xxx)
        source_decay = self._equation_term_source_decay(xxx, ttt)
        for i in range(len(self.c_source)):
            y_term = self._equation_term_y(i, xxx, yyy)
            cxyt_step = 1 / 8 * self.c_source[i] * source_decay * (x_term + additional_x) * y_term * z_term
            cxyt += cxyt_step
        if self._mode == "instant_reaction":
            self.cxyt_noBC = cxyt.copy()
            cxyt -= self.biodegradation_capacity
            cxyt = np.where(cxyt < 0, 0, cxyt)
        self.has_run = True
        return cxyt


# bioscreen ; Domenico - additional term


class Bioscreen(Anatrans):
    """Model class using the analytical solution implemented in the BIOSCREEN screening model, Newell et al. (1997).

    This model is an exact implementation of the transport equations implemented in the BIOSCREEN screening model of
    Newell et al. (1997), which is based on the Domenico (1987) analytical model. Using a truncated version of the
    equation used in the Anatrans model. This model is implemented as a method of comparison with the original BIOSCREEN
    software. And is included for legacy reasons, since it is the first model implemented in the mibitrans package,
    serving as a basis for the other models. However, caution should be taken when using this model, since a varying
    error is introduced by using the truncated analytical solution. The error is most prominent for shorter times and
    distances from the source, and depends on the ratio of flow velocity and longitudinal dispersivity. For modelling,
    the Anatrans (untruncated approximate solution) and Mibitrans (exact analytical solution) models are recommended
    instead.

    Domenico, P. A. (1987). An analytical model for multidimensional transport of a decaying contaminant species.
    Journal of Hydrology, 91(1-2), 49-58.

    Newell, C. J., McLeod, R. K., & Gonzales, J. R. (1997). BIOSCREEN natural attenuation decision support
    system version 1.4 revisions, Tech. rep., U.S. EPA.
    """

    def __init__(
        self,
        hydrological_parameters,
        attenuation_parameters,
        source_parameters,
        model_parameters,
        verbose=False,
    ):
        """Initialize model object.

        Args:
            hydrological_parameters (mibitrans.data.parameters.HydrologicalParameters) : Dataclass object containing
                hydrological parameters from HydrologicalParameters.
            attenuation_parameters (mibitrans.data.read.AttenuationParameters) : Dataclass object containing adsorption,
                degradation and diffusion parameters from AttenuationParameters.
            source_parameters (mibitrans.data.read.SourceParameters) : Dataclass object containing source parameters
                from SourceParameters.
            model_parameters (mibitrans.data.read.ModelParameters) : Dataclass object containing model parameters from
                ModelParameters.
            verbose (bool, optional): Verbose mode. Defaults to False.

        Attributes:
            mode (str) : Current model mode. Is 'linear' by default. Once instant reaction parameters are provided. Use
                this attribute to switch between 'linear' and 'instant_reaction' models.
            cxyt (np.ndarray) : Output array containing concentrations in model domain, in [g/m^3]. Indexed as [t,y,x]
            relative_cxyt (np.ndarray) : Output array with concentrations in model domain, divided by the maximum source
                zone concentration at t=0. Indexed as [t,y,x].
            x (np.ndarray) : Discretized model x-dimension, in [m].
            y (np.ndarray) : Discretized model y-dimension, in [y].
            t (np.ndarray) : Discretized model t-dimension, in [days].
            c_source (np.ndarray) : Nett source zone concentrations, accounting for source superposition, in [g/m^3].
            vr (float) : Retarded groundwater flow velocity, in [m/d].
            k_source (float) : Source zone decay rate, in [1/days].
            biodegradation_capacity (float) : Maximum capacity of biodegradation given provided electron acceptor
            concentrations, in [g/m^3].

        Methods:
            run : Run model with current parameters, returns Results object.
            sample : Calculate concentration at any given position and point in time.
            instant_reaction : Activate the instant reaction model by providing electron acceptor concentrations. And
                optionally electron acceptor utilization factors. Switch between model modes by using the mode
                attribute.

        Raises:
            TypeError : If input is not of the correct Dataclass.

        """
        super().__init__(hydrological_parameters, attenuation_parameters, source_parameters, model_parameters, verbose)

    @property
    def short_description(self):
        """Short description of model type."""
        return "Bioscreen model"

    def _calculate_concentration_for_all_xyt(self, xxx, yyy, ttt):
        # Difference with the Anatrans solution is the lack of additional term.
        cxyt = 0
        with np.errstate(divide="ignore", invalid="ignore"):
            decay_sqrt = np.sqrt(1 + 4 * self._decay_rate * self._hyd_pars.alpha_x / self.rv)
            x_term = self._equation_term_x(xxx, ttt, decay_sqrt)
            z_term = self._equation_term_z(xxx)
            source_decay = self._equation_term_source_decay(xxx, ttt)
            for i in range(len(self.c_source)):
                y_term = self._equation_term_y(i, xxx, yyy)
                cxyt_step = 1 / 8 * self.c_source[i] * source_decay * x_term * y_term * z_term
                cxyt += cxyt_step
        if self._mode == "instant_reaction":
            self.cxyt_noBC = cxyt.copy()
            cxyt -= self.biodegradation_capacity
            cxyt = np.where(cxyt < 0, 0, cxyt)
        self.has_run = True
        return cxyt
