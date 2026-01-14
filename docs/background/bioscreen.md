# BIOSCREEN

## General
[*BIOSCREEN*](https://www.epa.gov/water-research/bioscreen-natural-attenuation-decision-support-system) has been developed by the U.S. Environmental Protection Agency (EPA) in collaboration with the U.S. Air Force as natural attenuation decision support tool. Its is meant as screening tool to determine if a full-scale evaluation of a contaminated site is needed. Thereby it is not replacing a more involved numerical model, but serves as preliminary step to evaluate the necessity of an involved numerical model. 

*BIOSCREEN* calculates contaminant concentration distributions in 3D for a constant source under uniform flow conditions based on the advection dispersion equation. The source domain can contain several zones of different concentrations. Three modes of decay can be chose to represent degradation processes:
* no decay
* linear decay
* instantaneous biodegradation reaction

Input parameters for *BIOSCREEN* are few compared to numerical models. After parameter entry, visualization of concentrations at the plume centreline or as a 3D plume can be performed for each of the three decay modes. Version 1.4 provides the option to calculate mass balances of the plume and the source. *BIOSCREEN* comes with two sets of example data.

*BIOSCREEN* is implemented in Excel, and provides a graphical interface. The latest version is *BIOSCREEN* 1.4. It was released to be compatible with Microsoft Excel 5.0. Analysis is performed using macro scripts and cellular calculations. Calculations inside the Excel sheets are hidden in the background, and give insight into the actual calculations behind the model.
The spatial resolution of the plume is fixed to 11 steps over the plume length and 5 steps lateral to the plume. Temporal resolution is fixed to 10 time steps. This resolution is relative to the set model extent and time.

## Transport model

Transport is modeled based on the three-dimensional advection-dispersion equation (ADE) with linear equilibrium adsorption for uniform flow in $x$-direction:

$$
\begin{equation}\tag{1}\label{eq:01_ADE3D}
    R\frac{\partial C}{\partial t} = -v\frac{\partial C}{\partial x} + D_{x}\frac{\partial ^2 C}{\partial x^2} + D_y \frac{\partial^2 C}{\partial y^2} + D_z \frac{\partial^2 C}{\partial z^2} + r_{sinks}
\end{equation}
$$

Here $C(x,y,z,t)$ is the contaminant concentration in space $(x,y,z)$ and time $t$, $R$ is the linear equilibrium retardation factor, $v$ is the uniform groundwater velocity in $x$-direction, $D_{x}$, $D_{y}$ and $D_{y}$ are the longitudinal, transverse horizontal and transverse vertical dispersion coefficients, respectively. $r_{sinks}$ represents the sink term as result to decay/degradation. $r_{sinks} = 0$ represents no decay/degradation. For linear decay, $r_{sinks} = -\lambda C$. In the equation $\frac{\partial C}{\partial t}$ represents the change of the concentration over time. $-v\frac{\partial C}{\partial x}$ is the change of concentration in the direction of the groundwater gradient due to advection. $D_{x}\frac{\partial ^2 C}{\partial x^2} + D_{y}\frac{\partial ^2 C}{\partial y^2}+ + D_{z}\frac{\partial^2 C}{\partial z^2}$ represent the change in concentration to dispersion. 

A specific solution of the ADE $\eqref{eq:01_ADE3D}$, depends on the specific form of the sink-term $r_{sinks}$ as well as initial and boundary conditions.
Transport in *BIOSCREEN* is modelled based on the analytical model of *Domenico*, [1987]. 

### Domenico Model

An analytical model for the ADE for initial and boundary conditions representing typical spill field situation has been presented by *Domenico, 1987*. The solution has not been derived mathematically rigorous as analytical solution of the ADE, but has been composed of analytical solutions for the individual processes. As *Domenico*, [1987] writes, *an exact solution to this problem cannot avaid some form of numerical integration*. The provided analytical expression approximates the concentration distribution of a decaying species that is released to the aquifer as an extended pulse. *West et al., 2007* provides a detailed overview on the effects of the approximations in the *Domenico*-solution and potential error it can introduce to solute transport predictions.

Specifically, the *Domenico*-model takes the following assumption on initial and boundary conditions:
* Flow is uniform in $x$ direction with constant velocity $v$.
* The contaminant decays continuously at a rate of $\lambda$ (independent of position).
* There is no adsorption/retardation of the contaminant.
* The contaminant is released to the aquifer within a source plane of width $Y$, height $Z$, located at $x=0$ and centered at $y=0$ and $z=0$, so the center of plume is always located along the $x$-axis.
* The input of contaminant at the source is constant over time with an amount of $C_0$ that does not change over time. Specifically, it assumes that the source is not subject to depletion or internal decay/degradation that reduced source concentrations.

The equation describing the solute distribution in time and space of the *Domenico*-model reads:

$$
\begin{align}\tag{2}
    C(x, y, z, t) &= \frac{C_{0}}{8} \exp \left[ \frac{x\left(1-\sqrt{1+4\lambda \alpha_x/v}\right)}{2\alpha_x}\right] \\
    &\quad \cdot \operatorname{erfc} \left[ \frac{x - vt\sqrt{1+4\lambda \alpha_x/v}}{2\sqrt{\alpha_x vt }} \right] \\
    &\quad \cdot \left\{ \operatorname{erf} \left[ \frac{y + Y/2}{2\sqrt{\alpha_y x}} \right] - \operatorname{erf} \left[ \frac{y - Y/2}{2\sqrt{\alpha_y x)}} \right] \right\} \\
    &\quad \cdot \Biggl. \left\{ \operatorname{erf} \left[ \frac{z +Z/2}{2\sqrt{\alpha_z x)}} \right] - \operatorname{erf} \left[ \frac{z-Z/2}{2\sqrt{\alpha_z x}} \right] \right\} \nonumber
\end{align}
$$

where $C(x,y,z,t)$ is the contaminant concentration in $M/V^3$ at position $(x,y,z)$ and time $t$. $v$ is the groundwater flow velocity in $m/d$. $\alpha_x$ is the longitudinal dispersivity (in the x-direction). $\alpha_y$ is the transverse horizontal dispersivity (in the y-direction). $Y$ and $Z$ are the width and thickness/height of the source in the saturated zone in $m$.

The first two terms account for the transport of contaminants in horizontal direction due to advection and longitudinal dispersion while being subject to continuous decay. This part is based on the plug flow model of *Bear, 1979* and is in line with other derivations of analytical solutions for the 1D ADE assuming uniform flow, constant decay and continuous input of contaminant.

The two terms in the second row account for the dilution of the plume due to transverse dispersion. They were derived by *Domenico and Palciauskas, 1982*. The model is formulated as a boundary value problem that approximates the spreading. However, the terms are not in line with mathematically rigorous derived analytical solutions for this process. The terms $\sqrt{\alpha_{y/z} x)}$ in the denominator should read $\sqrt{\alpha_{y/z} vt)}$. Consequently, this form effectively assumes a linear increase of lateral dispersion as function of the plume travel distance. This is a misinterpretation of the concept of transverse spreading. Lateral spreading takes place at every distance. Note, that this should not be mixed up with the characteristics of longitudinal dispersion that is evolving over distance up to a asymptotic value that depends on aquifer heterogeneity (i.e. also not linearly with plume travel distance).

Bioscreen makes use of the *Domenico*-model in various adaptions (*Newell et al., 1996, 1997*).


### BIOSCREEN Implementation of Decay and No-decay solutions

Bioscreen does not directly use the Domenico model, but an extension regarding source handling and retardation. For the cases of linear decay of the contaminant ($\lambda \neq 0$) and the case of no decay ($\lambda = 0$), they include a source decay term accounting for reduction/depletion of the input concentration from the source in the form of $C_0(x,t) = C_0 \exp{\left( -k_s \left( t-\frac{x}{v}\right)\right)}$.

The analytical expression used by BIOSCREEN for no-decay ($\lambda = 0$) and linear decay ($\lambda \neq 0$) reads (*Newell et al., 1997*):

$$
\begin{align}\tag{3}
    C(x, y, z, t) & = \frac{C_{0}}{8}\exp{\left( -k_s \left( t-\frac{x}{v}\right)\right)} \\ 
    &\quad \cdot  \exp \left[ \frac{x\left(1-\sqrt{1+4\lambda \alpha_x/v}\right)}{2\alpha_x}\right] 
    \cdot \operatorname{erfc} \left[ \frac{x - vt\sqrt{1+4\lambda \alpha_x/v}}{2\sqrt{\alpha_x vt }} \right] \\
    &\quad \cdot \left\{ \operatorname{erf} \left[ \frac{y + Y/2}{2\sqrt{\alpha_y x}} \right] - \operatorname{erf} \left[ \frac{y - Y/2}{2\sqrt{\alpha_y x)}} \right] \right\} \\
    &\quad \cdot \Biggl. \left\{ \operatorname{erf} \left[ \frac{Z/2}{2\sqrt{\alpha_z x)}} \right] - \operatorname{erf} \left[ \frac{Z/2}{2\sqrt{\alpha_z x}} \right] \right\} \nonumber
\end{align}
$$

Note that here $z=0$ compared to Eq. $\eqref{eq:02_domenico}$, indicating that Bioscreen only evaluates the solution $z=0$ and does not provide a vertically resolved solution.

Bioscreen also uses the Domenico model in Eq. $\eqref{eq:02_domenico}$ as starting point for developing the model for the instantaneous reaction. 

$$
\begin{align}\tag{4}
    C(x, y, z, t) &= \frac{C_{0}}{8}\exp{\left( -k_s \left( t-\frac{x}{v}\right)+BC\right)} \\
    &\quad \cdot \operatorname{erfc} \left[ \frac{x - vt}{2\sqrt{\alpha_x vt }} \right] \\
    &\quad \cdot \left\{ \operatorname{erf} \left[ \frac{y + Y/2}{2\sqrt{\alpha_y x}} \right] - \operatorname{erf} \left[ \frac{y - Y/2}{2\sqrt{\alpha_y x)}} \right] \right\} \\
    &\quad \cdot \Biggl. \left\{ \operatorname{erf} \left[ \frac{Z/2}{2\sqrt{\alpha_z x)}} \right] - \operatorname{erf} \left[ \frac{Z/2}{2\sqrt{\alpha_z x}} \right] \right\} -BC \nonumber
\end{align}
$$

#### BIOSCREEN Implementation of Instant reaction model

BIOSCREEN uses a superposition approach in combination with the *Domenico-model* to model instantaneous aerobic and anaerobic reactions in groundwater, based on available concentrations of electron acceptors (EAs). They argue that the comparison to more complex models resolving the processes shows good agreement (within the range of assumptions) justifying the use of this heuristic approach for reactive transport modelling.

#### Principle

The general model idea is:
* Calculate how much contaminant can be consumed by an instantaneous reactions with all present EA’s based on their global concentrations. That amount is called biodegredation capacity BC. The BC is a lumped value that reflects the potential contaminant mass removal of available EAs.
* Calculate the spatially distributed contaminant concentration as if there is no decay. But adapt the handling of the source concentration: the source zone concentration is the sum of the measured source zone concentration and the biodegradation capacity BC.
* Subtract the BC from the calculated concentration for every location and time, i.e. reduce the concentration by what can be consumed by reactions.

This procedure is implemented as a superposition of the reaction to the Domenico model. By this method, contaminant mass concentrations are transported
conservatively and then corrected at any location and time within the flow field by subtracting 1 mg/L organic mass for each mg/L of BC provided by all of the available electron acceptors.

Specifically, the equation for the instant reaction model in BIOSCREEN reads:

$$
\begin{align}\tag{5}
    C(x, y, t) + BC &= \sum_{i=0}^{n} \Biggl\{ \left( C^*_{0,i} \exp \left[-k_s^{inst} \left(t - \frac{xR}{v} \right)\right] + BC \right) \bigr. \\
    &\quad \quad \quad \cdot \left\{ \frac{1}{8} \operatorname{erfc} \left[ \frac{x - \frac{vt}{R}}{2\sqrt{\alpha_x \frac{vt}{R}}} \right] \right\}  \\
    &\quad \quad \quad \cdot \left\{ \operatorname{erf} \left[ \frac{y + Y^*_i}{2\sqrt{\alpha_y x}} \right] - \operatorname{erf} \left[ \frac{y - Y^*_i}{2\sqrt{\alpha_y x)}} \right] \right\} \\
    &\quad \quad \quad \cdot \Biggl. \left\{ \operatorname{erf} \left[ \frac{Z}{2\sqrt{\alpha_z x)}} \right] - \operatorname{erf} \left[ \frac{-Z}{2\sqrt{\alpha_z x}} \right] \right\} \Biggr\} 
\end{align}
$$

#### Biodegradation Capacity BC

The Biodegradation capacity is calculated via:

$$
\begin{equation}\tag{6}
   BC = \sum_{O,N,S} (\bar C_i^\mathrm{upgradient} - C_i^\mathrm{source})/UF_i + \sum_{Fe2+,CH_4^+}  \bar C_j^\mathrm{source}/UF_j 
\end{equation}
$$

here $\bar C_i^\mathrm{upgradient}$ is the average upgradient concentrations and $C_i^\mathrm{source})$ is the minimum source concentration of $i=$  oxygen, nitrate, sulfate. $\bar C_j^\mathrm{source}$ is the average source concentration of $j=Fe2+,CH_4^+$. 
$UF_i$ and $UF_j$ are the utilization factors for each EA that was developed based on the stoichimetric ratios of the reactions (see below).

The usage of measured concentrations is based on the following assumption for the different EAs:
* Available EA from oxygen, nitrate and sulfate are transported with the groundwater, thus their concentrations are replenished.Their available concentrations for biodegradation are equal to the difference between average upgradient and minimum source concentrations (which is considered the available background concentration). A key assumptions here is that GW upstream is unaffected and their is full consumption of these AE in source zone.
* Available EA from iron-reducing and methanogenesis reactions are determined from measure concentrations of metabolic by-products: ferrous iron (Fe2+) and methane in source zone. EAs of these reaction are difficult to quantify (C02 is produced as end product of other reactions and ferric iron (Fe3+) is dissolved from the aquifer matrix. A key assumption here is that the model does not account for depletion of $Fe3+$ in aquifer matrix.

The calculated value of BC for the given global concentrations of EAs/EA-byproducts provides an estimate of the biodegradation capacity of the groundwater flowing through the source zone and plume and the aquifer soil matrix. In the instantaneous reaction model, it is assumed that EA's are consumed to full capacity for contaminant reduction. This includes the assumption that all reactions occur over the entire area of the contaminant plume.

#### Utilization factors UF

Utilization factors are based on the stoichimetric ratios of the reactions between the four BTEX compounds and each of the EAs as presented in *Wiedemeier, 1995*. The utilization factor for oxygen, nitrate, and sulfate can be developed showing the stoichiometric ratio of EA consumed to the mass of dissolved hydrocarbon degraded in the biodegradation reactions. Utilization factors for iron reduction and methanogenesis can be developed from the ratio of generated mass of metabolic by-products to mass of dissolved hydrocarbon degraded. 


| **EA/Byproduct**   | **UF (gm/gm)** |
|--------------------|----------------|
| Oxygen             | 3.14           |
| Nitrate            | 4.9            |
| Sulfate            | 4.7            |
| Ferrous Iron       | 21.8           |
| Methane            | 0.78           |

*Table: BTEX utilization factors (UF) for redox reactions*

Note that UFs are limited to reactions of EAs with BTEX constituents. When aiming to model other contaminants, the utilization factors would need to be adapted. Alternatively, available oxygen, nitrate, iron, sulfate, and methane concentrations could be adjusted accordingly to reflect alternate utilization factors. 

### References

Bear, J., Hydraulics of groundwater, London ; New York : McGraw-Hill International Book Co., 1979

[Domenico, P. A., and V. V. Palciauskas, Alternative Boundaries in Solid Waste Management, Groundwater, 20 (3), 303–311, 1982](https://doi.org/10.1111/j.1745-6584.1982.tb01351.x)

[Domenico, P., An analytical model for multidimensional transport of a decaying contaminant species, Journal of Hydrology, 91, 49–58, doi:10.1016/0022-1694(87)90127-2, 1987.] (https://doi.org/10.1016/0022-1694(87)90127-2)

[Newell, C. J., R. K. Mcleod, J. R. Gonzales, and J. T. Wilson, BIOSCREEN natural attenuation decision support system user’s manual version 1.3, Tech. rep., U.S. EPA, 1996.](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P1007K50.TXT)

Newell, C. J., R. K. Mcleod, J. R. Gonzales, and J. T. Wilson, BIOSCREEN Natural Attenuation Decision Support System Version 1.4 Revisions, U.S. EPA, 1997.

[West, M. R., B. H. Kueper, and M. J. Ungs, On the use and error of approximation in the Domenico (1987) solution, Groundwater, 45 (2), 126–135, 2007](https://doi.org/10.1111/j.1745-6584.2006.00280.x)

Wiedemeier, T. H., J. T. Wilson, D. H. Kampbell, R. N. Miller, and J. E. Hansen, Technical protocol for implementing intrinsic remediation with long-term monitoring for natural attenuation of fuel contamination dissolved in groundwater. Volume II, Tech. Rep. AD-A–324247/6/XAB, Parsons Engineering Science, Inc., Denver, CO (United States), 1995

## Bugs in *BIOSCREEN*

During the setup of `mibitrans`, *BIOSCREEN* has been thoroughly tested. It was found to be host to a minor erroneous calculation for the instant reaction model. *BIOSCREEN* uses calculations from the *no decay model* and corrects them for the different source decay coefficient and source zone concentrations. However, in these corrections, the wrong source decay coefficient is used, resulting in an underestimation of modelled biodegradation. The size of the error is determined by choice of parameters relating to source decay and biodegradation capacity. 

