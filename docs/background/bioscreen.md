# BIOSCREEN

## General

[*BIOSCREEN*](https://www.epa.gov/water-research/bioscreen-natural-attenuation-decision-support-system) has been developed by the U.S. Environmental Protection Agency (EPA) in collaboration with the U.S. Air Force as natural attenuation decision support tool *Newell et al, (1996, 1997)*. It is meant as screening tool to determine if a full-scale evaluation of a contaminated site is needed. Thereby it is not replacing a more involved numerical model, but serves as preliminary step to evaluate the necessity of an involved numerical model. 

*BIOSCREEN* calculates contaminant concentration distributions in 3D for a constant source under uniform flow conditions using an empirical solution of the advection dispersion equation. The source domain can contain several zones of different concentrations. Three modes of decay can be chose to represent degradation processes:

* no decay
* linear decay
* instantaneous biodegradation reaction

Version 1.4 provides the option to calculate mass balances of the plume and the source. *BIOSCREEN* comes with two sets of example data.

*BIOSCREEN* is implemented in Excel, and provides a graphical interface. The latest version is *BIOSCREEN* 1.4. It was released to be compatible with Microsoft Excel 5.0. Analysis is performed using macro scripts and cellular calculations. Calculations inside the Excel sheets are hidden in the background, and give insight into the actual calculations behind the model.
The spatial resolution of the plume is fixed to 11 steps over the plume length and 5 steps lateral to the plume. Temporal resolution is fixed to 10 time steps. This resolution is relative to the set model extent and time. *BIOSCREEN* only evaluates the solution z=0 and does not provide a vertically resolved solution.

---

## Transport Model

Transport in *BIOSCREEN* is modelled based on the analytical model of *Domenico*, [1987]. The actual implementation is the same as the equation for the *Bioscreen* Model Class in the `mibitrans` package. 
The *Domenico*-solution has not been derived with mathematical rigor as analytical solution of the ADE, but has been composed of analytical solutions for the individual processes. As *Domenico*, [1987] writes, *"an exact solution to this problem cannot avoid some form of numerical integration"*. The provided analytical expression approximates the concentration distribution of a decaying species that is released to the aquifer as an extended pulse. *West et al., 2007* provides a detailed overview on the effects of the approximations in the *Domenico*-solution and potential error it can introduce to solute transport predictions.

Specifically, the *Domenico*-model takes the following assumption on initial and boundary conditions:

* Flow is uniform in $x$ direction with constant velocity $v$.
* The contaminant decays continuously at a rate of $\lambda$ (independent of position).
* There is no adsorption/retardation of the contaminant.
* The contaminant is released to the aquifer within a source plane of width $Y$, height $Z$, located at $x=0$ and centered at $y=0$ and $z=0$, so the center of plume is always located along the $x$-axis.
* The input of contaminant at the source is constant over time with an amount of $C_0$ that does not change over time. Specifically, it assumes that the source is not subject to depletion or internal decay/degradation that reduced source concentrations.

*BIOSCREEN* does not directly use the Domenico model, but an extension regarding source handling and retardation. For the cases of linear decay of the contaminant (and the case of no decay), a source decay term is included that accounts for reduction/depletion of the input concentration from the source in the form.

*BIOSCREEN* uses a superposition approach in combination with the *Domenico-model* to model instantaneous aerobic and anaerobic reactions in groundwater, based on available concentrations of electron acceptors (EAs). They argue that the comparison to more complex models resolving the processes shows good agreement (within the range of assumptions) justifying the use of this heuristic approach for reactive transport modelling. This *Instant Reaction Model* approach can be understood as a special form of the fringe-degradation approach. 

---

## Instantaneous Reaction Model

### Background / History

The *instantaneous reaction model* is based on the concept developed by *Connor et al., [1994]* using a solution for contaminant transport (advection, dispersion, retardation) and adapting it to account for the influence of chemical reactions during transport. Specifically, the effect of degradation is superimposed based on domain average electron acceptor (EA) concentrations. *Connor et al., [1994]* compared his method with results of the BIOPLUME II model *[Rifai et al., 1990]* and determined that results are similar for biodegradable contaminants with retardation factors $R \leq 6$.

---

### Concept

The *instantaneous reaction model* is based on the assumption that degradation is fast, i.e. instantaneous compared to the time scale of other processes, such as transport. Specifically, the rate of biodegradation is not limited by microbial kinetics. For oxygen this is well accepted and shown by *Borden et al., [1986]*. For anaerobic degradation this is concluded from results of *Connor et al., [1994]*. The model further assumes that contaminants and EA's travel at the same transport rate, i.e. retardation is negligible. The validity of this assumption was shown by *Borden et al., [1986]* for oxygen and some hydrocarbon, but not for other EA's. *[Newell et al., 1996]* reports that the assumption holds for other EAs and retardation factors  $R\leq6$.

Starting point of the model is calculating the quantity of contaminant that can be consumed through biodegradation with all present EA's based on their global concentrations. This quantity is called the *biodegradation capacity* $BC$. Details on the calculation of $BC$ are given below. $BC$ is a lumped value that reflects the potential contaminant mass removal of available EAs through all biodegradation reactions, including aerobic and anaerobic reactions during transport.

In a second step, the spatial concentration distribution is calculated by combining $BC$ with an analytical model of the solute transport equation (*mibitrans*, *anatrans*, *bioscreen*) for no decay.
When

$$
\begin{equation}\tag{1}
C(x, y, z, t) = \frac{C_s}{8} f(x,y,z,t)
\end{equation}
$$

is an analytical solution for retarded, but non-decaying solute transport ($\mu = 0$), then the instantaneous reaction model is constructed via

$$
\begin{equation}\tag{2}
C_\mathrm{inst}(x, y, z, t)
=
\frac{(C_s + BC)}{8} f(x,y,z,t) - BC
\end{equation}
$$

Note that $C_s$ implies that source depletion as an option. For the bioscreen model, this reads e.g. $C_s = (C_0 + BC)\exp{\left( -\gamma_s \left( t-\frac{xR}{v}\right)\right)}$.

Eq. 2 for $C_\mathrm{inst}(x, y, z, t)$ shows that the source conditions are adapted by adding the biodegradation capacity $BC$. This reflects the source concentrations without any biodegradation. Consequently, the resulting transported and spatially distributed contaminant concentration is reduced with $BC$ over the entire domain. This reflects that calculated concentrations are reduced at every location and time by the capacity of biological activity to degrade contaminants, which is assumed to be constant.

---

### Biodegradation Capacity (BC)

The $BC$ is calculated from EA concentrations by

$$
BC =
\sum_{j=O_2,NO_3^-,SO_4^{2-}}
\frac{\bar C_j^\mathrm{upgradient} - C_j^\mathrm{source}}{UF_j}
+\sum_{k=Fe^{2+},CH_4}
\frac{\bar C_j^\mathrm{source}}{UF_k}
$$

here $\bar C_j^\mathrm{upgradient}$ is the average upgradient concentrations and $C_j^\mathrm{source}$ is the minimum source concentration of electron acceptors $j=O_2, NO_3^-, SO_4^{2-}$. $\bar C_k^\mathrm{source}$ is the average source concentration of the metabolic by-products $k=Fe^{2+},CH_4$. 
$UF_j$ and $UF_k$ are the utilization factors for each EA that were developed by *Wiedemeier et al., [1995]* based on the stoichiometric ratios of the reactions (see below).

| EA / Byproduct | UF ($gm/gm$) |
| -------------- |--------------|
| Oxygen         | 3.14         |
| Nitrate        | 4.9          |
| Sulfate        | 4.7          |
| Ferrous Iron   | 21.8         |
| Methane        | 0.78         |

BTEX utilization factors (UF) for redox reactions.

---

### Utilization Factors (UF)

Utilization factors as derived by *Wiedemeier et al., [1995]* and used by *Bioscreen* with values given in previous Table are based on the degradation reactions, here exemplified for benzene:

$$
C_6H_6 + 7.5O_2 \to 6CO_2 + 3H_2O
$$

$$
6NO_3^- + 6H^+ + C_6H_6 \to 6CO_2 + 6H_2O + 3N_2
$$

$$
60H^+ + 30Fe(OH)_3 + C_6H_6 \to 6CO_2 + 6Fe^{2+} + 78H_2O
$$

$$
7.5H^+ + 3.75SO_4^{2-} + C_6H_6 \to 6CO_2 + 3.75H_2S + 3H_2O
$$

$$
C_6H_6 + 4.5H_2O \to 2.25CO_2 + 3.75CH_4
$$

The utilization factor for oxygen, nitrate, and sulfate can be developed showing the stoichiometric ratio of EA consumed to the mass of dissolved hydrocarbon degraded in the biodegradation reactions. Utilization factors for iron reduction and methanogenesis can be developed from the ratio of generated mass of metabolic by-products to mass of dissolved hydrocarbon degraded. Values represented in the Table are averages of stoichiometric ratios for all four BTEX constituents.

Note that UFs are limited to reactions of EAs with BTEX constituents. When aiming to model other contaminants, the utilization factors would need to be adapted. Alternatively, available oxygen, nitrate, iron, sulfate, and methane concentrations could be adjusted accordingly to reflect alternate utilization factors. 

---

### Transport Model Implementation

*BIOSCREEN* calculates contaminant concentrations for the *instantaneous reaction model* by modifying the analytical model of *Domenico*, [1987] (same equation as the `Bioscreen` Model Class in the `mibitrans`) without decay ($\mu = 0$, thus $P = 1$) as follows:

$$
\begin{aligned}
C(x, y, t)
&=
\sum_{i=1}^{n}
\left\{
\left[
\frac{C^*_{0,i} + BC_n}{8}
\exp\left(
-\gamma_s \left( t - \frac{x}{v_R} \right)
\right)
\right]
\operatorname{erfc}
\left(
\frac{x - v_R t}{2\sqrt{\alpha_x v_R t}}
\right)
\right.
\\
&\quad \cdot
\left[
\operatorname{erf}
\left(
\frac{y + Y_i}{2\sqrt{\alpha_y x}}
\right) -
\operatorname{erf}
\left(
\frac{y - Y_i}{2\sqrt{\alpha_y x}}
\right)
\right]
\\
&\quad \cdot
\left.
\left[
\operatorname{erf}
\left(
\frac{Z}{2\sqrt{\alpha_z x}}
\right)-
\operatorname{erf}
\left(
\frac{-Z}{2\sqrt{\alpha_z x}}
\right)
\right]
\right\}
-BC
\end{aligned}
$$

Here, $BC_n$ indicates that for a model with multiple source zones, the biodegradation capacity is only added to the outermost source zone. These adaptation are the same for the *anatrans* and *mibitrans* equations.

---

### Instantaneous Reaction with Source Depletion

The instant reaction model assumes that biodegradation takes place inside the source zone. Thus, the source depletion coefficient is calculated differently as

$$
\gamma_s^{\mathrm{inst}}
=
\frac{Q \cdot (\bar C_0 + BC)}{m_{\mathrm{s,0}}}
$$

Here, the biodegradation capacity $BC$ is added to the average source zone contaminant concentration, resulting in a greater contaminant flux out of the source zone, and thus, faster source depletion.

This has also implication on the calculation of the mass balance. For the model option *instant reaction*, the current mass of pure phase contaminant in the source zone at any point in time $t$, is calculated via

$$
m_{\mathrm{s}}(t) = m_\mathrm{s,0}\cdot e^{-\gamma_{s}^{inst} \cdot t}
$$

---


## Minor Errors in *BIOSCREEN*

During the setup of `mibitrans`, *BIOSCREEN* has been thoroughly tested. It was found to be host to a minor erroneous calculation for the instant reaction model. *BIOSCREEN* uses calculations from the *no decay model* and corrects them for the different source decay coefficient and source zone concentrations. However, in these corrections, the wrong source decay coefficient is used, resulting in an underestimation of modelled biodegradation. The size of the error is determined by choice of parameters relating to source decay and biodegradation capacity.


---

## References

Borden, R. C., and P. B. Bedient, *Transport of dissolved hydrocarbons influenced by oxygen-limited biodegradation: 1. Theoretical development*, Water Resources Research, 22 (13), 1973–1982, doi:10.1029/WR022i013p01973, 1986.

Connor, J., C. Newell, J. Nevin, and H. Rifai, *Guidelines for use of groundwater spreadsheet models in risk-based corrective action design*, Tech. rep., National Ground Water Association, Proceedings of the Petroleum Hydrocarbons and Organic Chemicals in Ground Water Conference, Houston, Texas, 1994.

[Domenico, P., An analytical model for multidimensional transport of a decaying contaminant species, Journal of Hydrology, 91, 49–58, doi:10.1016/0022-1694(87)90127-2, 1987.] (https://doi.org/10.1016/0022-1694(87)90127-2)

[Newell, C. J., R. K. Mcleod, J. R. Gonzales, and J. T. Wilson, BIOSCREEN natural attenuation decision support system user’s manual version 1.3, Tech. rep., U.S. EPA, 1996.](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P1007K50.TXT)

[Newell, C. J., R. K. McLeod, and J. R. Gonzales, BIOSCREEN natural attenuation decision support system version 1.4 revisions, Tech. rep., U.S. EPA, 1997.](https://d3pcsg2wjq9izr.cloudfront.net/files/6377/download/651291/0-3.pdf)

Rifai, H. S., and P. B. Bedient, *Comparison of biodegradation kinetics with an instantaneous reaction model for groundwater*, Water Resources Research, 26 (4), 637–645, doi:10.1029/WR026i004p00637, 1990

[West, M. R., B. H. Kueper, and M. J. Ungs, On the use and error of approximation in the Domenico (1987) solution, Groundwater, 45 (2), 126–135, 2007](https://doi.org/10.1111/j.1745-6584.2006.00280.x)

Wiedemeier, T. H., J. T. Wilson, D. H. Kampbell, R. N. Miller, and J. E. Hansen, *Technical protocol for implementing intrinsic remediation with long-term monitoring for natural attenuation of fuel contamination dissolved in groundwater*. Volume II, Tech. Rep. AD-A–324247/6/XAB, Parsons Engineering Science, Inc., Denver, CO (United States), 1995.
