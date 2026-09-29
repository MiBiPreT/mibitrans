## Core-fringe degradation model

The concept of a combined core and fringe analytical degradation model is first introduced in Gutierrez-Neri et al. (2009).
Distinction is made between processes happening at the plume fringes, and in the plume interior (core).


### Background
Borden et al. (1986), uses a numerical model of the 2D-ADE for a non-degrading, non-reactive tracer to model transport of both electron acceptor (EA) and electron donor (ED). 
Provided velocity field and transport coefficients are the same, a superposition method can be used to determine resulting concentration distribution of a degradation reaction between EA and ED. 
Defined as $C_{ED} = C_{ED}^T - \frac{C_{EA}^T}{f}$.
Where $C_{ED}$ is the electron donor concentration after biodegradation, 
$C_{ED}^T$ is the total (pre-biodegradation) electron donor concentration, 
$C_{EA}^T$ is the total (pre-biodegradation) electron acceptor concentration 
and f is the mass ratio of EA to ED consumed in the biodegradation reaction.
This equation implies that $1mol$ of EA is responsible for degradation of $1mol$ ED.
This is however not the case in most biodegradation reactions, in which case concentrations of ED and EA should be presented in the same stoichiometric units, which is elaborated upon in `Biodegradation stoichiometry` below.
For simplicity going forward, assume a conceptual biodegradation reaction $n_1 ED + n_2 EA \to n_3 P$, where $n$ is the stoichiometric coefficient and P is the degradation product, where $n$ is 1.

In case of Borden et al. (1986) specifically, Oxygen was used as electron acceptor, but this method of superposition is only considerd valid if the reaction between the EA and ED is rapid compared to the groundwater flow velocity, where EA concentration is the only limiting factor.
Under this assumption, the biodegradation reaction is considered as an 'instantaneous' reaction, where EA and ED cannot co-exist.
Furthermore, for multiple EA species, it is assumed that there is no preference in the order in which they occur.
And that there is no preference for a specific pair of ED and EA.

### Analytical solution fringe degradation

This approach was included in the BIOPLUME software, and later as the instantaneous reaction model in BIOSCREEN, where the same principle was applied to an analytical solution.
The paper of Koussis et al. (2003) investigates the validity of the instantaneous reaction model, by comparing it to a kinetic numerical model.
Concluding that in general it is valid except near the source zone or during initial plume development.
Though, no analytical derivation of this approach was published until the paper of Ham et al. (2004).
Discussing the mathematical derivation is somewhat beyond the scope, but importantly, the analysis of Ham et al. (2004) also highlights further boundary values and limitations of the approach.
For the derivation of the solution, continuous injection at the origin $(x, y) \to (0,0)$ is assumed and $t \to \infty$.
Due to the latter boundary condition, the derived analytical equation is for steady state, and for determining the plume length.
Later in the paper, considerations are discussed for a non-steady state plume. For transient cases, the approach only valid where the retardation factor of the ED is equal to that of the EA.
Though, as the retardation factor only affects transient development of contaminant plume, this assumption is not required for the steady state solution (Ham et al., 2004).

### Combined core-fringe degradation
In the case that reaction rate is considered the limiting factor instead of the electron acceptor concentration, the assumption of an instant reaction no longer holds.
Specifically for anaerobic degradation like interactions with metal oxides or fermentation, this can be the case, usually in the plume center, where oxygen is largely absent.
These processes can then be approximated by first-order, linear decay.
[_Claimed in this manner by Gutierrez-Neri et al. (2009), which cites Lovley et al., (1989), however, no support of this is present in that paper. Though, other literature does take this approach for anaerobic degradation as well_].
Assumption for modelling the biodegradation as first-order decay is that reaction rate is indeed dependent on electron donor concentration, and that the electron acceptor does not deplete.
The paper of Gutierrez-Neri et al. (2009) describes an analytical solution for the combination of core degradation and fringe degradation.
In general, this solution looks like:
$$ \tag{1}
C_{ED}(x,y,z,t) = F_{ED}(x,y,z,t, C^0_{ED},\lambda) - (C^0_{EA} - F_{EA}(x,y,z,t,C^0_{EA}))
$$
where $F_{ED}$ is a transport equation for a linearly decaying species with source concentration of $C^0_{ED}$ and degradation rate $\lambda$, 
and $F_{EA}$ is that same transport equation for a conservative species with source concentration $C^0_{EA}$.
All other transport parameters are the same in $F_{ED}$ and $F_{EA}$.
For $F$, Gutierrez-Neri et al. (2009) uses the solution of Domenico and Robbins (1985), with linear degradations as in Domenico (1987).
Hunkeler et al. (2010) pointed out a few issues in the overall formulation of the core-fringe degradation equation, where eq.7 reads:

$$ \tag{2}
C_{ED}(x,y,t) = 
\left\{ 
    \begin{array}{l}
        0 \text{ for } C^T_{ED}\cdot \left(K(x,\lambda) \cdot\frac{F_1(x,\lambda,t)}{F_1(x,t)} + \frac{C^0_{EA}}{C^0_{ED}}\right) \leq C_{0,EA}&\\
         C^T_{ED}\cdot \left(K(x,\lambda) \cdot\frac{F_1(x,\lambda,t)}{F_1(x,t)} + \frac{C^0_{EA}}{C^0_{ED}}\right) - C_{0,EA} \text{ elsewhere}
    \end{array}
\right\}
$$
Where

$$\tag{3}
C^T_{ED} = \frac{C^0_{ED}}{4}\cdot F_1(x,t) \cdot F_2(x,y)
$$

$$\tag{4}
K(x,\lambda) = \exp\left[ \left( \frac{x}{2\alpha_x} \right) \left( 1-\sqrt{1+\frac{4\lambda \alpha_x}{v}} \right) \right]
$$

$$\tag{5}
F_1(x,t) = \text{erfc}\left(\frac{x-vt}{2\sqrt{\alpha_x-vt}} \right)
$$

$$\tag{6}
F_1(x,\lambda,t) = \text{erfc}\left(\frac{x-vt\sqrt{(1+4\alpha_x/v)}}{2\sqrt{\alpha_x-vt}} \right)
$$

$$\tag{7}
F_2(x,y) = \left[ \text{erf}\left(\frac{y + Y/2}{2 \sqrt{\alpha_yx}}\right) - \text{erf}\left( \frac{y - Y/2}{2\sqrt{\alpha_yx}} \right) \right]
$$

### Biodegradation stoichiometry
As is, equation 1 is only valid if 1g of EA is used to degrade 1g of ED, which is not the case in any degradation reaction.
For any equation relating EA and ED concentrations to be valid, there has to be unit stoichiometry.
In case of BIOSCREEN, this is encapsuled in the Biodegradation Capacity (BC). 
Here electron acceptor concentrations are divided by their respective utilization factor, which is the mass ratio in the biodegradation reaction between the specific EA and ED. i.e. grams of EA required to degrade 1g of ED.
As BIOSCREEN is specifically designed for modelling BTEX degradation, the mass ratio for an EA is the average between B, T, E and X. 
By summing the values for each EA, one gets the BC, which expresses EA concentrations in the form of concentration of degradable ED.
The approach taken by Gutierrez-Neri et al. (2009) instead considers $C_{ED}$ and $C_{EA}$ in mol electrons per unit of volume.
This then also causes unit stoichiometry, as for each electron donated, an electron is accepted.
Conversion to mol electrons per volume is achieved by considering the half-reactions of the ED and EA. 

### Implementation to _mibitrans_
While Gutierrez-Neri et al. (2009) and Hunkeler et al. (2010) use the Domenico solution for the core-fringe model, the governing principles apply to the exact (Mibitrans) solution as well.
However, due to the integral involved in the exact solution, it does not resolve to the form of equation 2, like the fully analytical solutions would.
Instead, _mibitrans_ uses the implemented transport equations as is, and calculates ED and EA concentrations as follows:
$$ \tag{8}
C_{ED}(x,y,t)  = \operatorname{max} \left[F_{ED}(x,y,t) -  \left(C_{0,EA} - F_{EA}(x,y,t)\right),0 \right] \\
$$
$$ \tag{9}
C_{EA}(x,y,t)  = \operatorname{max} \left[  \left(C_{0,EA} - F_{EA}(x,y,t) \right) - F_{ED}(x,y,t),0 \right]
$$

Here, $F$ is the `Mibitrans`, `Anatrans` or `Bioscreen` model equation, as described in the `Model implementations` section of the documentation.
$F_{ED}$ uses the source concentration of the ED and the input decay rate.
$F_{EA}$ uses the source concentration of the EA and does not decay.

Four additional parameters are needed for this implementation of fringe degradation; electron acceptor concentrations, electron acceptor molar weights, electron acceptor stoichiometric ratio and electron donor molar weight.
Where electron acceptor stoichiometric ratio is referring to stoichiometry of EA:ED of the entire (not half) biodegradation reaction.
Parameters for any number electron acceptors can be given, which are transformed based on the stoichiometric parameters given and subsequently summed to obtain a single $C_{0,EA}$.

### (Source) superposition
The paper of Gutierrez-Neri et al. (2009) mentions different source conditions (source depletion, pulse-injection), but does not involve multiple source zones (source-superposition), as implemented in _mibitrans_.
As both core-fringe degradation and multiple source zones are based on the principle of superposition, and are linear.
Therefore, it can be argued that (at least in the mathematical sense), it is valid to combine these methods.
However, no mathematical evidence to support this has been found in literature.

### Sources

Borden, R. C., Bedient, P. B., Lee, M. D., Ward, C. H., & Wilson, J. T. (1986). Transport of dissolved hydrocarbons influenced by oxygen-limited biodegradation: 2. Field application. Water Resources Research, 22(13), 1983–1990. https://doi.org/10.1029/WR022i013p01983

Gutierrez-Neri, M., Ham, P. A. S., Schotting, R. J., & Lerner, D. N. (2009). Analytical modelling of fringe and core biodegradation in groundwater plumes. Journal of Contaminant Hydrology, 107(1), 1–9. https://doi.org/10.1016/j.jconhyd.2009.02.007

Ham, P. A. S., Schotting, R. J., Prommer, H., & Davis, G. B. (2004). Effects of hydrodynamic dispersion on plume lengths for instantaneous bimolecular reactions. Advances in Water Resources, 27(8), 803–813. https://doi.org/10.1016/j.advwatres.2004.05.008

Hunkeler, D., Höhener, P., & Atteia, O. (2010). Comments on “Analytical modelling of fringe and core biodegradation in groundwater plumes.” by Gutierrez-Neri et al. in J. Contam. Hydrol. 107: 1–9. Journal of Contaminant Hydrology, 117(1), 1–6. https://doi.org/10.1016/j.jconhyd.2010.06.009

Koussis, A. D., Pesmajoglou, S., & Syriopoulou, D. (2003). Modelling biodegradation of hydrocarbons in aquifers: When is the use of the instantaneous reaction approximation justified? Journal of Contaminant Hydrology, 60(3), 287–305. https://doi.org/10.1016/S0169-7722(02)00083-9
