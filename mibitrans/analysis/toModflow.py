import os
import tempfile
import flopy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


class MibitransToModflow:
    """Set up and run a MODFLOW model based on Mibitrans model parameters.

    This class utilizes the Dataclass objects from a Mibitrans model to define
    parameters and build and run the MODFLOW model (Langevin et al., 2026).

    Langevin, C. D., Hughes, J. D., Provost, A. M., Russcher, M. J., Niswonger, R. G., Panday, S., Merrick, D., Morway,
    E. D., Reno, M. J., Bonelli, W. P., Boyce, S. E., & Banta, E. R. (2026). MODFLOW 6 modular hydrologic model
    (Version 6.8.0.dev0) [Computer software]. U.S. Geological Survey. https://doi.org/10.5066/F76Q1VQV
    """

    def __init__(
        self,
        hydrological_parameters,
        attenuation_parameters,
        source_parameters,
        model_parameters
    ):
        """Initialize model object.

        Args:
            hydrological_parameters (mibitrans.data.parameters.HydrologicalParameters) : Dataclass object containing
                hydrological parameters, dispersion and diffusion from HydrologicalParameters.
            attenuation_parameters (mibitrans.data.parameters.AttenuationParameters) : Dataclass object containing
                adsorption and degradation parameters from AttenuationParameters.
            source_parameters (mibitrans.data.parameters.SourceParameters) : Dataclass object containing source
                parameters from SourceParameters.
            model_parameters (mibitrans.data.parameters.ModelParameters) : Dataclass object containing model parameters
                from ModelParameters.
        """
        self.hydro = hydrological_parameters
        self.att = attenuation_parameters
        self.source = source_parameters
        self.model = model_parameters

        # Set up the file names
        self.exe_name_mf = "mf2005"
        self.exe_name_mt = "mt3dms"

        self.mf = None
        self.mt = None

        self.conc = None
        self.times = None

        self._parameters()

    def _parameters(self):
        """Calculate and set model parameters.

        Calculates hydrological, attenuation, model, source and advection parameters from the model input data
        and stores them as attributes of the model object.

        Hydrological parameters:
            v (float) : Groundwater velocity [m/d].
            hk (float) : Hydraulic conductivity [m/d].
            prsity (float) : Effective porosity of porous medium [-].
            al (float) : Longitudinal dispersivity [m].
            trpt (float) : Ratio of horizontal transverse dispersivity to longitudinal dispersivity [-].
            trpv (float) : Ratio of vertical transverse dispersivity to longitudinal dispersivity [-].
            i (float) : Hydraulic gradient [m/m].
            laytyp (int) : Layer type; either confined or convertible [-].

        Attenuation parameters:
            lambda1 (float) : First-order reaction rate [d^-1].
            retardation (float) : Retardation factor [-].
            rhob (float) : Bulk density of aquifer medium [g/m3].
            kd (float) : Distribution coefficient [m3/g].
            isothm (int) : Represents the type of sorption used in the model [-].
            ireact (int) : Represents the type of kinetic rate reaction used in the model [-].
            igetsc (int) : Represents whether the initial concentration of the nonequilibrium sorbed or the
                immobile phase of all the species should be read [-].

        Source parameters:
            icbund (int) : Array containing the boundary condition for each cell [-].
            sconc (float) : Array containing the starting concentration for each grid cell [g/m3].

        Model parameters:
            perlen (float) : Array containing length of the current stress period [d].
            nstp (int) : Number of time steps within each stress period [-].
            delr (float) : Array containing the interval lengths along the rows [-].
            delc (float) : Array containing the interval lengths along the columns [-].
            nrow (int) : Total number of rows [-].
            ncol (int) : Total number of columns [-].
            nlay (int) : Total number of layers [-].
            top (float) : Array containing the top elevation of the first layer [m].
            botm (float) : Array containing the bottom elevation for each cell [m].
            ipakcb (int) : Represents whether cell-by-cell budget data should be saved [m].
            ibound (int) : Array containing a value representing constant head, inactive cell,
                or variable head for each grid cell [-].
            strt (float) : Array containing the value of the starting head for each grid cell [m].

        Advection parameters:
            dceps (float) : Represents value below which advective transport is considered negligible (RCCG) [-].
            nplane (int) : Variable which indicates whether a random of fixed pattern is selected for the
                initial placement of moving particles [-].
            npl (int) : Number of initial particles placed in cells in which RCCG is less or equal to dceps [-].
            nph (int) : Number of initial particles placed in cells in which RCCG is greater than dceps [-].
            npmin (int) : Maximal amount of particles present in a cell [-].
            npmax (int) : Minimal amount of particles present in a cell [-].
            mixelm (int) : Integer variable to distinguish between different advection solution options [-].
            percel (float) : Number of cells, or the fraction of a cell in which advection is allowed in any
                direction in one transport step [-].
        """
        # Hydrological parameters
        self.v = self.hydro.velocity # [m/d]
        self.hk = self._get_parameter("hydraulic conductivity", self.hydro.h_conductivity, "m/d") # [m/d]
        self.prsity = self.hydro.porosity # [-]
        self.al = self.hydro.alpha_x # [m]
        self.trpt = self.hydro.alpha_y/self.al # [-]
        self.trpv = self.hydro.alpha_z/self.al # [-]

        if self.hydro.h_gradient is not None:
            self.i = self.hydro.h_gradient
        else:
            self.i = (self.v * self.prsity)/self.hk # [m/m]

        self.laytyp = 0 # Confined layer [-]

        # Attenuation parameters
        self.lambda1 = self.att.decay_rate # [1/d]
        self.rhob = self._get_parameter("bulk density", self.att.bulk_density, "g/m3") # [g/m3]

            # Retardation [m]
        if self.att.retardation is not None:
            self.retardation = self.att.retardation
        elif self.att.partition_coefficient is not None and self.att.fraction_organic_carbon is not None:
            self.retardation = (
                1
                + self.rhob
                * self.att.partition_coefficient
                * self.att.fraction_organic_carbon
                / self.prsity
            )

            # Distribution coefficient [m3/g]
        if self.att.partition_coefficient is not None and self.att.fraction_organic_carbon is not None:
            self.kd = self.att.partition_coefficient * self.att.fraction_organic_carbon
        else:
            self.kd = (self.retardation - 1.0) * self.prsity / self.rhob

        self.isothm = 1 # Linear isotherm [-]
        self.ireact = 1 # First-order irreversible reaction [-]
        self.igetsc = 0 # Initial concentration of the sorbed or immobile phase is not read [-]

        # Model parameters
        delv = self.source.depth # [m]
        self.perlen = self.model.model_time # [d]
        self.nstp = self.perlen/self.model.dt # [-]
        self.delr = self.model.dx # [m]
        self.delc = self.model.dy # [m]
        self.nrow = int((2*self.model.model_width)/self.delc) # [-]
        self.ncol = int(self.model.model_length/self.delr) # [-]
        self.nlay = 1 # [-]
        self.top = delv # [m]
        self.botm = [0] # [m]
        self.ipakcb = 53 # Cell-by-cell budget data is saved [-]

        self.ibound = np.ones((self.nlay, self.nrow, self.ncol), dtype=int) # [-]
        self.ibound[0, :, 0] = -1 # Constant head [-]
        self.ibound[0, :, -1] = -1 # Constant head [-]
        self.strt = np.zeros((self.nlay, self.nrow, self.ncol), dtype=float) # [m]
        lx = (self.ncol-1) * self.delr # [m]
        h1 = self.i * lx # [m]
        self.strt[0, :, 0] = h1 + delv # [m]
        self.strt[0, :, -1] = delv # [m]

        # Source parameters
        sbound = self.source.source_zone_boundary # [m]
        szoneconc = self.source.source_zone_concentration # [g/m3]
        amount_zones = len(sbound) # [-]
        self.icbund = np.ones((self.nlay, self.nrow, self.ncol), dtype=int) # [-]
        self.icbund[0, :, 0] = -1 # Constant concentration [-]
        self.sconc = np.zeros((self.nlay, self.nrow, self.ncol), dtype=float) # [g/m3]
        y = (np.arange(self.nrow) + 0.5)*self.delc # [m]
        y_center = (self.nrow*self.delc)/2 # [m]
            # Determining which cells fall within each source zone and assigning the concentration
        for i in range(amount_zones - 1, -1, -1):
            mask = np.abs(y - y_center) <= sbound[i] # boolean [-]
            self.sconc[0,mask,0] = szoneconc[i] # [g/m3]

        # Advection parameters
        self.dceps = 1.0*10**-5 # [-]
        self.nplane = 2 # 3D-simulations [-]
        self.npl = 0 # [-]
        self.nph = 4 # [-]
        self.npmin = 0 # [-]
        self.npmax = 8 # [-]
        self.mixelm = 2 # Backward-tracking modified method of characteristics (MMOC) [-]
        self.percel = 0.5 # [-]

        self._write_report()

    def _write_report(self):
        """Writes a short report after parameter initialization including assumptions and calculated values."""
        assumptions = []

        if self.hydro.alpha_z == 1e-10:
            assumptions.append("Transverse vertical dispersivity = 1e-10 m")
        if self.lambda1 == 0:
            assumptions.append("Decay rate = 0 [1/d]")
        if self.retardation == 1:
            assumptions.append("Retardation = 1")

        calculations = []

        if self.hydro.h_gradient is None:
            calculations.append(f"Hydraulic gradient = {self.i} m/m")
        calculations.append(f"Distribution coefficient = {self.kd} m3/g")
        if self.att.partition_coefficient is not None and self.att.fraction_organic_carbon is not None:
            calculations.append(f"Retardation = {self.retardation}")
        calculations.append(f"Hydraulic head left = {self.strt[0, 0, 0]} m")
        calculations.append(f"Hydraulic head right = {self.strt[0, 0, -1]} m")

        report = ("Initialization completed.")

        if assumptions:
            report += "\nThe following assumptions have been made:\n"
            report += "\n".join(f"   - {assumption}" for assumption in assumptions)

        report += "\nThe following parameters have been calculated:\n"
        report += "\n".join(f"   - {calculation}" for calculation in calculations)

        print(report)

    def _get_parameter(self, name, value, unit):
        """Check if parameter exists in the given input data classes.
        
        Checks whether the parameter is part of the input data classes and prompts the user to
        provide the value if it is not.

        Args:
            name (str) : The name of the requested parameter.
            value (float) : The value of the requested parameter.
            unit (str) : The unit corresponding to the requested parameter.

        Returns:
            value (float): The value corresponding to the requested parameter.
        """
        if value is None:
            print(f"Enter {name} in {unit}, using a period as the decimal separator")
            value = float(input())
        return value

    def to_modflow(self, model_ws = None, temporary = True):
        """Generate MODFLOW and MT3DMS model objects.

        Args:
            model_ws (str, optional) : Path to the working directory of the MODFLOW model.
                If not specified, a temporary working directory is used.
            temporary (bool, optional) : Indicates whether a temporary working directory is used for the model.
        """
        if model_ws is not None:
            self.model_ws = model_ws
        elif temporary:
            self._temp_dir = tempfile.TemporaryDirectory(prefix="mt3dms_")
            self.model_ws = self._temp_dir.name
        else:
            self.model_ws = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "model"
            )
        os.makedirs(self.model_ws, exist_ok=True)

        self._generate_modflow()
        self._generate_mt3dms()

        print(f"Model objects have been saved to {self.model_ws}")

    def _generate_modflow(self):
        """Builds a MODFLOW model object.

        Packages:
            - DIS (Discretization Package Clas)
            - BAS (Basic Package Class)
            - LPF (Layer Property Flow Package Class)
            - PCG (Preconditioned Conjugate-Gradient Package class)
            - LMT (Link-MT3DMS Package Class)

        Returns:
            mf (model) : MODFLOW model object.
        """
        # Defining the name of the model and making the MODFLOW model
        modelname_mf = "modflow_mf"
        self.mf = flopy.modflow.Modflow(
            modelname = modelname_mf,
            model_ws = self.model_ws,
            exe_name = self.exe_name_mf,
            )

        # Adding all necessary MODFLOW Packages to the model
        flopy.modflow.ModflowDis(
            self.mf,
            nlay=self.nlay,
            nrow=self.nrow,
            ncol=self.ncol,
            delr=self.delr,
            delc=self.delc,
            top=self.top,
            botm=self.botm,
            perlen=self.perlen,
            nstp = self.nstp,
        )
        flopy.modflow.ModflowBas(self.mf, ibound=self.ibound, strt=self.strt)
        flopy.modflow.ModflowLpf(self.mf, ipakcb=self.ipakcb, hk=self.hk, laytyp=self.laytyp)
        flopy.modflow.ModflowPcg(self.mf)
        flopy.modflow.ModflowLmt(self.mf)

    def _generate_mt3dms(self):
        """Builds a MT3DMS model object.

        Packages:
            - BTN (Basic Transport Package Class)
            - ADV (Advection Package Class)
            - DSP (Dispersion Package Class)
            - RCT (Chemical Reaction Package Class)
            - SSM (Source And Sink Mixing Package Class)
            - GCG (Generalized Conjugate Gradient Package Class)

        Returns:
            mt (model) : MT3DMS model object.
        """
        # Defining the name and making the MT3DMS model (with the MODFLOW model as part of the input)
        modelname_mt = "mt3dms_mt"
        self.mt = flopy.mt3d.Mt3dms(
            modelname=modelname_mt,
            model_ws=self.model_ws,
            exe_name=self.exe_name_mt,
            modflowmodel=self.mf,
        )

        # Adding all necessary MT3DMS Packages to the model
        flopy.mt3d.Mt3dBtn(self.mt, icbund=self.icbund, prsity=self.prsity, sconc=self.sconc)
        flopy.mt3d.Mt3dAdv(
            self.mt,
            mixelm=self.mixelm,
            dceps=self.dceps,
            nplane=self.nplane,
            npl=self.npl,
            nph=self.nph,
            npmin=self.npmin,
            npmax=self.npmax,
            nlsink=self.nplane,
            npsink=self.nph,
            percel=self.percel,
        )
        flopy.mt3d.Mt3dDsp(self.mt, al=self.al, trpt=self.trpt, trpv=self.trpv)
        flopy.mt3d.Mt3dRct(
            self.mt,
            isothm=self.isothm,
            ireact=self.ireact,
            igetsc=self.igetsc,
            rhob=self.rhob,
            sp1=self.kd,
            rc1=self.lambda1,
            rc2=self.lambda1,
        )
        flopy.mt3d.Mt3dSsm(self.mt)
        flopy.mt3d.Mt3dGcg(self.mt)

    def run_modflow(self, verbose=True):
        """Calculate the concentration for all discretized z, y, x and t using MODFLOW and MT3DMS.

        Returns:
            conc (float) : Array containing the concentration for all grid cells and simulation times [g/m3]
            times (float) : Array containing all of the simulation times [d]
        """
        self.mf.write_input()
        self.mf.run_model(silent=True)
        if verbose:
            print("The Modflow model run is completed.")

        self.mt.write_input()
        fname = os.path.join(self.model_ws, "MT3D001.UCN")
        if os.path.isfile(fname):
            os.remove(fname)
        self.mt.run_model(silent=True)
        if verbose:
            print("The MT3DMS model run is completed.")

        if verbose:
            print("The results are being retrieved.")
        fname = os.path.join(self.model_ws, "MT3D001.UCN")
        ucnobj = flopy.utils.UcnFile(fname)
        self.times = ucnobj.get_times()
        self.conc = ucnobj.get_alldata()
        if verbose:
            print("The results have been retrieved succesfully.")

        return self.conc, self.times

    def centerline_modflow(self, label=None):
        """Plot center of contaminant plume at the end of the simulation.

        Args:
            label (str, optional) : Label used to identify the data series in the graph legend.
        """
        # Setting the size of figure
        mpl.rcParams["figure.figsize"] = (8, 8)

        row = self.mf.modelgrid.nrow // 2
        x = self.mf.modelgrid.xcellcenters[row,:]
        y = self.conc[-1, 0, row, :]
        t = self.perlen

        plt.plot(x, y, label = label)
        plt.xlabel("Distance from source [m]")
        plt.ylabel("Concentration[g/m3]")
        plt.title(f"Centerline plot of MODFLOW model, at t = {t} days")
        if label is not None:
            plt.legend()
