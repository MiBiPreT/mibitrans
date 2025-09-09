# Theoretical Background on Subsurface Transport Models

## Analytical Models: No Decay & Linear Decay

### Transport Equation

Transport is modelled based on the two-dimensional **advection-dispersion equation** (ADE) with linear equilibrium adsorption: 
\begin{equation}
    R\frac{\partial C}{\partial t} = -v\frac{\partial C}{\partial x} + D_{x}\frac{\partial ^2 C}{\partial x^2} + D_y \frac{\partial^2 C}{\partial y^2} + r_{sinks}
\end{equation}

Here $C(x,y,t)$ is the contaminant concentration in space $(x,y)$ and time $t$, $R$ is the linear equilibrium retardation factor, $v$ is the uniform groundwater velocity in $x$-direction, $D_{x}$ and $D_{y}$ are the longitudinal and transverse horizontal dispersion coefficients, respectively. $r_{sinks}$ represents the sink term as result to decay/degredation. Given the specific choice of degredation model, $r_{sinks}$ takes different mathematical expression, including $r_{sinks} = 0$ for no decay.

In the equation $\frac{\partial C}{\partial t}$ represens the change of the concentration over time. $-v\frac{\partial C}{\partial x}$ is the change of concentration in the direction of the groundwater gradient due to advection. $D_{x}\frac{\partial ^2 C}{\partial x^2} + D_{y}\frac{\partial ^2 C}{\partial y^2}$ represent the change in concentration to dispersion. 

### Degradation Models in Transport Equation

Three cases of degradation during solute transport are considered: 
* no decay ($r_{sinks} = 0$)
* linear decay $r_{sinks} = - \mu C$ or $r_{sinks} = - \mu R C$
* instantaneous biodegradation reaction

For linear decay, the sink terms reads $r_{sinks} = \mu R C$, where $\mu$ is the first order decay constant. $R$ is included when the adsorbed contaminant still undergoes biodegradation.

For the case of decay as instantaneous biodegradation reaction, there is no explicit representation of $ $r_{sinks}$, but the analytical solution of the ADE is modified differently. The instantaneous biodegradation reaction model will be discussed separately.

### Initial and Boundary Conditions

We consider contaminant transport with a continues source as typical for LNAPL and DNAPL contaminations. The source is modelled as plane source located at the center of the coordinate system $x = 0$. It can consist of multiple zones perpendicular to the groundwater flow, i.e. along the $y$-axis. 
The corresponding mathematical specifications for the ADE are given by the initial conditions
\begin{equation}
    C(x=0,Y_i,t = 0) = C_{0,i} \quad C(x!=0,y,z,t = 0) = 0
\end{equation}
where $C_{0,i}$ is the initial concentration within region $i$ along the $y$-axis. Everywhere else, the contaminant concentration is initially zero.

The transport domain is considered infinite. The corresponding boundary condition reads:
\begin{equation}
    lim_{r\to \infty} C(r,t)= 0
\end{equation}
where $r = \sqrt{x^2 + y^2}$ is the infinite radial distance at which concentration is zero.

The contaminant source is seen as a separate, pure phase of contaminant dissolving into the groundwater that also undergoes decay. Over time, the pure phase diminishes due to transport of contaminated groundwater out of the source zone. The processes is represented as an exponential decay of the source zone concentrations. Here, $k_s$ is the source decay coefficient, which is dependant on the contaminant flux out of the source zone. The corresponding boundary condition reads:
\begin{equation}
    C_{source,i}(t)= C_{0,i} \exp \left(-k_s \left( t-\frac{x}{v} \right) \right)
\end{equation}

In the case of *no decay*, the set of initial and boundary conditions as given above is complete. 
The set of initial and boundary conditions is also representing the case of *linear decay*. In this case, it is assumed that the modelled biodegradation does not take place inside the source zone. Biodegradation only starts after leaving the source zone and has no effects on the contaminant concentrations in the source zone self. So, source decay refers here to the source depletion due to dissolution not a decay within the source as result of biodegradation

The *instant reaction model* is not a direct solution of the ADE with the given initial and boundary conditions. It is formed from an alteration of the analytical solution for the *no decay* model. Consequently, underlying assumptions are more involved.

### Assumptions
Note that the transport model represented by the ADE and its initial and boundary conditions is based on several assumtions:
* Groundwater flow is uniform solely in $x$-direction (i.e. the coordinate system is aligned with the major flow direction).
* Dispersion is following Scheideggers Law. Groundwater flow velocity is fast enough for molecular diffusion can be ignored ($D_{disp} >> D_{eff}$).
* The aquifer is homogeneous and isotropic.
* The contaminant source is given in the form of a line source, so dilution effects in transverse vertical direction can be neglected (and the transport situation is simplyfied to 2D.
* Adsorption is represented by a linear isotherm and reversible.
* The source zone is symmetrical in $x = 0$, which highest concentrations in the central zone around $y = 0$.

### Analyical solutions

An analytical solution for the given ADE under the specified initial and boundary conditions has been presented by *Domenico*, [1986].

#### No decay

The equation implemented in `mibitrans` for the *no decay* model reads:

\begin{equation}
    C(x, y, t) &= \sum_{i=0}^{n} \Biggl\{ \left( C^*_{0,i} \exp \left[-k_s \left(t - \frac{xR}{v} \right)\right] \right) \bigr. \\
    &\quad \quad \quad \cdot \left\{ \frac{1}{8} \operatorname{erfc} \left[ \frac{x - \frac{vt}{R}}{2\sqrt{\alpha_x \frac{vt}{R}}} \right] + \frac{1}{2} \exp \left[ \frac{xv}{\alpha_xR}\right] \cdot \operatorname{erfc}\left[ \frac{x + \frac{vt}{R}}{2\sqrt{\alpha_x \frac{vt}{R}}}\right]  \right\}  \\
    &\quad \quad \quad \cdot \left\{ \operatorname{erf} \left[ \frac{y + Y^*_i}{2\sqrt{\alpha_y x}} \right] - \operatorname{erf} \left[ \frac{y - Y^*_i}{2\sqrt{\alpha_y x)}} \right] \right\} \\
    &\quad \quad \quad \cdot \Biggl. \left\{ \operatorname{erf} \left[ \frac{Z}{2\sqrt{\alpha_z x)}} \right] - \operatorname{erf} \left[ \frac{-Z}{2\sqrt{\alpha_z x}} \right] \right\} \Biggr\} \\
\end{equation}
where $C(x,y,t)$ is the contaminant concentrationin $M/V^3$ at position $(x,y)$ and time $t$. $\sum_{i=0}^{n}$ is the sum of plume concentrations for the n amount of source zones i. $C^*_{o,i}$ is the nett initial concentration of contaminant at the source, in source zone i. $k_s$, is the first-order source decay coefficient in $days^{-1}$. $t$ is the time after contaminant release in $days$. $x$ is the distance down gradient from the source in $m$. $R$ is the retardation factor due to adsorption. $v$ is the groundwater flow velocity in $m/d$. $\alpha_x$ is the longitudinal dispersivity (in the x-direction) in $m$. $y$ is the position lateral to the source in $m$. $Y^*_i$ is the width of the nett source zone i in $m$. $\alpha_y$ is the transverse horizontal dispersivity (in the y-direction) in $m$. Lastly, $Z$ is the source thickness in the saturated zone in $m$. $\alpha_z$ is the transverse vertical dispersivity (in the z-direction) in $m$.

The model calculates the resulting contaminant plume for each source zone individually and then superimposes the plumes, resulting in the contaminant plume for the entire source. Concentrations in the source are used as the net source zone concentrations. Starting with the outer zone, going inwards, source zone concentrations from zones at the outside are subtracted from the source zone concentration.

The source decay coefficient $k_s$ is calculated from the flow parameters and source concentrations as
\begin{equation}
    k_{s} = \frac{Q \cdot C_{0,avg}}{m_{source,t=0}}
\end{equation}
where $Q$ is the volumetric flow rate in $m^3/d$, and $m_{source,t=0}$ is the initial source mass in $g$. $C_{0,avg}$ is the average source zone concentration in $g/m^3$. 

The volumentric flow rate $Q$ is calculated as
\begin{equation}
    Q =  v \cdot \theta_e \cdot Y \cdot Z
\end{equation}
where $\theta_e$ is the effective porosity. 

The averagee initial concentration $C_{0,avg}$ is calculated as
\begin{equation}
    C_{0,avg} = \frac{\sum(C_{0,i} \cdot Y_{i})}{Y}
\end{equation}
where $C_{0,i}$ is the initial source concentration of source zone i in $g/m^3$ and $Y_i$ is the width of the same source zone in $m$. Y is the total source width in $m$. It should be noted that in the event the source mass is considered to be infinite, the source decay equation for $k_s$ resolves to zero, and no source decay is considered. 

#### Linear decay

The implemented analytical solution for the *linear decay* model in `mibitrans` is given by:
\begin{equation}
    C(x, y, t) &= \sum_{i=0}^{n} \Biggl\{ C^*_{o,i} \exp\left[ -k_s \left(t - \frac{xR}{v}\right) \right] \Biggr. \\
    &\quad \quad \quad \cdot \frac{1}{8} \exp\left[\frac{x}{2\alpha_x} \left(1 - \sqrt{1 + \frac{4 \mu \alpha_x R}{v}}\; \, \right)\right] \\
    &\quad \quad \quad \cdot \operatorname{erfc}\left[\frac{x - \frac{v t}{R} \sqrt{1 + \frac{4 \mu \alpha_x R}{v} }}{2\sqrt{\alpha_x \frac{v t}{R}}} \; \right] \\
    &\quad \quad \quad \cdot \left\{ \operatorname{erf} \left(\frac{y + Y^*_i}{2\sqrt{\alpha_y x}}\right) - \operatorname{erf} \left(\frac{y - Y^*_i}{2\sqrt{\alpha_y x}}\right) \right\} \\
    &\quad \quad \quad \Biggl. \cdot \left\{ \operatorname{erf} \left(\frac{Z}{2\sqrt{\alpha_z x}}\right) - \operatorname{erf} \left(\frac{-Z}{2\sqrt{\alpha_z x}}\right) \right\} \Biggr\} \\
\end{equation}

where $\mu$ is the first order decay coefficient in $days^{-1}$. 

## The instant reaction model

The *instant reaction model* is not a direct solution of the ADE with the given initial and boundary conditions. It is formed from an alteration of the analytical solution for the *no decay* model. Consequently, underlying assumptions are more involved. 

It is assumed that biodegradation is a function of available electron acceptors, using the stoichiometric relations from degradation reactions. The biodegradation reactions are assumed to occur much faster than the replenishment of electron acceptors by groundwater flow. Therefore, biodegradation is seen as an instantaneous reaction relative to the flow velocity. This approach is applied field-wide and therefore, no zonation of biodegradation reactions is considered. Opposed to the linear decay model, the instant reaction model does assume biodegradation taking place in the source zone. One assumption is that concentrations everywhere in the domain represent the concentrations after degradation. Without degradation, the concentrations would have been elevated by an amount determined by the biodegradation capacity (BC), calculated from electron acceptor concentrations.

The mathematical expression for contaminant transport under the instant reaction model implemented in `mibitrans` reads:

\begin{equation}
C(x, y, t) + BC &= \sum_{i=0}^{n} \Biggl\{ \left( C^*_{0,i} \exp \left[-k_s^{inst} \left(t - \frac{xR}{v} \right)\right] + BC \right) \bigr. \\
    &\quad \quad \quad \cdot \left\{ \frac{1}{8} \operatorname{erfc} \left[ \frac{x - \frac{vt}{R}}{2\sqrt{\alpha_x \frac{vt}{R}}} \right] + \frac{1}{2} \exp \left[ \frac{xv}{\alpha_xR}\right] \cdot \operatorname{erfc}\left[ \frac{x + \frac{vt}{R}}{2\sqrt{\alpha_x \frac{vt}{R}}}\right]  \right\}  \\
    &\quad \quad \quad \cdot \left\{ \operatorname{erf} \left[ \frac{y + Y^*_i}{2\sqrt{\alpha_y x}} \right] - \operatorname{erf} \left[ \frac{y - Y^*_i}{2\sqrt{\alpha_y x)}} \right] \right\} \\
    &\quad \quad \quad \cdot \Biggl. \left\{ \operatorname{erf} \left[ \frac{Z}{2\sqrt{\alpha_z x)}} \right] - \operatorname{erf} \left[ \frac{-Z}{2\sqrt{\alpha_z x}} \right] \right\} \Biggr\} \\
\end{equation}

where $BC$ is the biodegradation capacity in $g/m^3$. $BC$ is added to the outermost net source zone concentration, and subsequently subtracted from the entire plume after superimposition. Calculation of the $BC$ is done by considering the stoichiometric ratio for biodegradation of mixed BTEX constituents, and electron acceptor concentrations:
\begin{equation}
BC = \sum \frac{C_{ea,n}}{U_n}
\end{equation}
where $C_{ea,n}$ is the concentration of the nth electron acceptor/byproduct involved in biodegradation in $g/m^3$. $U_n$ is the utilization factor of an electron acceptor/byproduct, as weight of electron acceptor/byproduct that is consumed/generated per weight of biodegraded contaminant. The electron acceptors used by this model are oxygen ($O_2$), nitrate ($NO_3^-$), ferric iron ($Fe^{3+}$), sulfate ($SO_4^{2-}$, and carbon dioxide ($CO_2$). The amount of $O_2$, $NO_3^-$ and $SO_4^{2-}$ consumed for biodegradation is calculated as the difference between the minimum source zone concentration and the average background concentration upgradient of the source zone ($\Delta O_2$, $\Delta NO_3^-$ and $\Delta SO_4^{2-}$). For $Fe^{3+}$ and $CO_2$, average concentrations of $Fe^{2+}$ and $CH_4$ are used as proxies instead. These are the reduced forms of $Fe^{3+}$ and $CO_2$ respectively, and here referred to as byproducts. 

The utilization factor is dependent on the stoichiometry of the biodegradation reaction, and thus on the contaminant type. Here, the same values for utilization factor are used as in *BIOSCREEN*. Values are based on the combined degradation of BTEX as reported in Wiedemeier et al., [1995]:
|Electron acceptor or byproduct | Utilization factor [-] |
|--------------------------------------------------------|
|Oxygen 			| 3.14 			 |
|Nitrate 			| 4.9 			 |	
|Sulfate 			| 4.7 			 |
|Ferrous Iron 			| 21.8 			 |
|Methane 			| 0.78 			 |

The instant reaction model assumes that biodegradation takes place inside the source zone, it uses a different calculation for the source decay coefficient:
\begin{equation}
    k_{s}^{inst} = \frac{Q \cdot (C_{0,avg} + BC)}{m_{source,t=0}}
\end{equation}
Here, the biodegradation capacity is added to the average source zone concentration, resulting in a greater contaminant flux out of the source zone, and thus, faster source decay. 


### References

[Domenico, P., An analytical model for multidimensional transport of a decaying contaminant species, Journal of Hydrology, 91, 49–58, doi:10.1016/0022-1694(87)90127-2, 1987.] (https://doi.org/10.1016/0022-1694(87)90127-2)

[Newell, C. J., R. K. Mcleod, J. R. Gonzales, and J. T. Wilson, BIOSCREEN natural attenuation decision support system user’s manual version 1.3, Tech. rep., U.S. EPA, 1996.](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P1007K50.TXT)

[Newell, C. J., R. K. McLeod, and J. R. Gonzales, BIOSCREEN natural attenuation decision support system version 1.4 revisions, Tech. rep., U.S. EPA, 1997.](https://d3pcsg2wjq9izr.cloudfront.net/files/6377/download/651291/0-3.pdf)

[Wiedemeier, T., J. T. Wilson, D. H. Kampbell, R. N. Miller, and J. E. Hansen, Technical protocol for implementing intrinsic remediation with long-term monitoring for natural attenuation of fuel contamination dissolved in groundwater volume II, Tech. rep., Air Force Center for Environmental Excellence, 1995.] (https://www.osti.gov/biblio/512639)

