# Introduction

## General

Contaminated sites pose a risk to humans and the environment. Innovative cleaning technologies are needed to remediate these sites and remove contaminants such as petroleum hydrocarbons (PHC), cyanides and hexachlorocyclohexane (HCH).

Conventional methods of contaminated site remediation are often costly and upkeep intensive. Bioremediation is an alternative of particular interest, as it degrades contaminants on-site. Assessment of ongoing biodegradation is an important step to check the feasibility for bioremediation. Similarly, modeling the fate of contaminants is key for understanding the processes involved and predicting bioremediation in the field. 

**Estimating the extend of a contaminant plume** through modelling is highly valuable for designing monitoring and site remediation strategies. Predicting the fate of contaminant in the subsurface requires combining simulations on groundwater flow, contaminant transport and chemical reactions. Simulations allows making predictions on amounts, locations and time scales of biodegradation as well as measures of bioremediation. Complex numerical models can provide a detailed picture by taking site geometry, complex flow pattern and contaminant source distrubution into account. But they require sufficient data and detailed knowledge of site conditions and come at high computational cost and expert knowledge for setup. Simple models, based on analytical and semi-analytical solutions of the subsurface transport equation do not necessarily provide a realistic distribution of the contaminant, but they allow a quick estimate of plume travel distances, plume extend and mass balance. They are thus a good first transport screening option.

The **purpose of the `mibitrans` package** is to provide such as transport screening model based on hydrogeological field data for biodegredation and bioremediation. 

## MIBIREM

[MIBIREM - Innovative technological toolbox for bioremediation](https://www.mibirem.eu/) is a EU funded consortium project by 12 international partners all over Europe working together to develop an *Innovative technological toolbox for bioremediation*. The project will develop molecular methods for the monitoring, isolation, cultivation and subsequent deposition of whole microbiomes. The toolbox will also include the methodology for the improvement of specific microbiome functions, including evolution and enrichment. The performance of selected microbiomes will be tested under real field conditions. The `mibitrans` package is part of this toolbox.

## Bioremediation

Bioremediation uses living organisms (including bacteria) to digest and neutralize environmental contaminants. Like the microbiome in the gut, which supports the body in digesting food, microbiomes at contaminated sites can degrade organic contaminant in soil and groundwater.

Processes relevant for general biodegradation and bioremediation prediction are:

+ hydrogeological flow and transport: this includes groundwater flow driven by hydraulic gradients, advective transport of contaminant, diffusion and dispersion
+ transformation and phase transition processes: dissolution, volatilization, adsorption/retardation, decay
+ biochemical processes: chemical reaction and microbial degradation
+ microbiome evolution: spatial distribution and temporal development of bacteria actively degrading contaminants under various and/or changing environmental conditions.

Modeling all these processes at the same time, requires a high level of model detail, spatially resolved parameter information and knowledge on initial and boundary conditions. This is typically not feasible in the field. Thus, we follow the approach to select and combine most relevant processes and have modeling sub-modules (repositories within the MiBiPreT organization) which can be used for data analysis and predictive modeling of individual or combined processes. At the same time, modules are designed to allow for coupling of processes and (modeling) sub-modules at a advanced stage of tool development.

## Functionality

`mibitrans` is supposed to serve as hydrogeological transport screening model based on field data for biodegredation and bioremediation. In parts `mibitrans` reproduces - and extends- the functionality of the *Excel*-based screening tool [`BIOSCREEN`](https://www.epa.gov/water-research/bioscreen-natural-attenuation-decision-support-system) [Newell et al., 1996]. Thus, transport is modelled based on the 3D advection dispersion equation considering linear equilibrium adsorption and various options for biodegradation. The option of no decay is also available. `mibitrans` is validated by comparing example field data with results from `BIOSCREEN`. 

## Structure

The core elements and folders for users of `mibitrans` are:

* The folder `mibitrans` contains the main functionality split up into folders for:
    * `data`
    * `transport`
    * `analysis` 
    * `visualization`
* The folder `examples` contains example workflows in the form of Jupyter-Notebooks outlining application of functionality on example data.


## References

[Newell, C. J., R. K. Mcleod, J. R. Gonzales, and J. T. Wilson, BIOSCREEN natural attenuation decision support system user’s manual version 1.3, Tech. rep., U.S. EPA, 1996.](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P1007K50.TXT)

[Newell, C. J., R. K. McLeod, and J. R. Gonzales, BIOSCREEN natural attenuation decision support system version 1.4 revisions, Tech. rep., U.S. EPA, 1997.](https://d3pcsg2wjq9izr.cloudfront.net/files/6377/download/651291/0-3.pdf)
