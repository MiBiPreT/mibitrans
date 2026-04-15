# BIOSCREEN

## General

[*BIOSCREEN*](https://www.epa.gov/water-research/bioscreen-natural-attenuation-decision-support-system) has been developed by the U.S. Environmental Protection Agency (EPA) in collaboration with the U.S. Air Force as natural attenuation decision support tool [Newell et al, 1996, 1997]. Its is meant as screening tool to determine if a full-scale evaluation of a contaminated site is needed. Thereby it is not replacing a more involved numerical model, but serves as preliminary step to evaluate the necessity of an involved numerical model. 

*BIOSCREEN* calculates contaminant concentration distributions in 3D for a constant source under uniform flow conditions using an empirical solution of the advection dispersion equation. The source domain can contain several zones of different concentrations. Three modes of decay can be chose to represent degradation processes:

* no decay
* linear decay
* instantaneous biodegradation reaction

Version 1.4 provides the option to calculate mass balances of the plume and the source. *BIOSCREEN* comes with two sets of example data.

*BIOSCREEN* is implemented in Excel, and provides a graphical interface. The latest version is *BIOSCREEN* 1.4. It was released to be compatible with Microsoft Excel 5.0. Analysis is performed using macro scripts and cellular calculations. Calculations inside the Excel sheets are hidden in the background, and give insight into the actual calculations behind the model.
The spatial resolution of the plume is fixed to 11 steps over the plume length and 5 steps lateral to the plume. Temporal resolution is fixed to 10 time steps. This resolution is relative to the set model extent and time. *BIOSCREEN* only evaluates the solution z=0 and does not provide a vertically resolved solution.

## Transport model

Transport in *BIOSCREEN* is modelled based on the analytical model of *Domenico*, [1987]. The actual implementation is the same as the equation for the *Bioscreen* Model Class in the `mibitrans` package (see section **Model Implementations**).

The *Domenico*-solution has not been derived mathematically rigorous as analytical solution of the ADE, but has been composed of analytical solutions for the individual processes. As *Domenico*, [1987] writes, *"an exact solution to this problem cannot avaid some form of numerical integration"*. The provided analytical expression approximates the concentration distribution of a decaying species that is released to the aquifer as an extended pulse. *West et al., 2007* provides a detailed overview on the effects of the approximations in the *Domenico*-solution and potential error it can introduce to solute transport predictions.

Specifically, the *Domenico*-model takes the following assumption on initial and boundary conditions:

* Flow is uniform in $x$ direction with constant velocity $v$.
* The contaminant decays continuously at a rate of $\lambda$ (independent of position).
* There is no adsorption/retardation of the contaminant.
* The contaminant is released to the aquifer within a source plane of width $Y$, height $Z$, located at $x=0$ and centered at $y=0$ and $z=0$, so the center of plume is always located along the $x$-axis.
* The input of contaminant at the source is constant over time with an amount of $C_0$ that does not change over time. Specifically, it assumes that the source is not subject to depletion or internal decay/degradation that reduced source concentrations.

*BIOSCREEN* does not directly use the Domenico model, but an extension regarding source handling and retardation. For the cases of linear decay of the contaminant (and the case of no decay), a source decay term is included that accounts for reduction/depletion of the input concentration from the source in the form.

BIOSCREEN uses a superposition approach in combination with the *Domenico-model* to model instantaneous aerobic and anaerobic reactions in groundwater, based on available concentrations of electron acceptors (EAs). They argue that the comparison to more complex models resolving the processes shows good agreement (within the range of assumptions) justifying the use of this heuristic approach for reactive transport modelling. Details of this model approach are outlined in section **Instant Reaction Model**.

## Minor errors in *BIOSCREEN*

During the setup of `mibitrans`, *BIOSCREEN* has been thoroughly tested. It was found to be host to a minor erroneous calculation for the instant reaction model. *BIOSCREEN* uses calculations from the *no decay model* and corrects them for the different source decay coefficient and source zone concentrations. However, in these corrections, the wrong source decay coefficient is used, resulting in an underestimation of modelled biodegradation. The size of the error is determined by choice of parameters relating to source decay and biodegradation capacity.

## References

[Domenico, P., An analytical model for multidimensional transport of a decaying contaminant species, Journal of Hydrology, 91, 49–58, doi:10.1016/0022-1694(87)90127-2, 1987.] (https://doi.org/10.1016/0022-1694(87)90127-2)

[Newell, C. J., R. K. Mcleod, J. R. Gonzales, and J. T. Wilson, BIOSCREEN natural attenuation decision support system user’s manual version 1.3, Tech. rep., U.S. EPA, 1996.](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P1007K50.TXT)

Newell, C. J., R. K. Mcleod, J. R. Gonzales, and J. T. Wilson, BIOSCREEN Natural Attenuation Decision Support System Version 1.4 Revisions, U.S. EPA, 1997.

[West, M. R., B. H. Kueper, and M. J. Ungs, On the use and error of approximation in the Domenico (1987) solution, Groundwater, 45 (2), 126–135, 2007](https://doi.org/10.1111/j.1745-6584.2006.00280.x)
