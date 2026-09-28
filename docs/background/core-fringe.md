## Core-fringe degradation model

The concept of a combined core and fringe analytical degradation model is first introduced in Gutierrez-Neri et al. (2009).
Distinction is made between processes happening at the plume fringes, and in the plume interior (core).


### Background
Borden et al. (1986), uses a numerical model of the 2D-ADE for a non-degrading, non-reactive tracer to model transport of both electron acceptor (EA) and electron donor (ED). 
Provided velocity field and transport coefficients are the same, a superposition method can be used to determine resulting concentration distribution of a degradation reaction between EA and ED. 
Defined as $C_{ED} = C_{ED}^T - \frac{C_{EA}^T}{F}$.
Where $C_{ED}$ is the electron donor concentration after biodegradation, 
$C_{ED}^T$ is the total (pre-biodegradation) electron donor concentration, 
$C_{EA}^T$ is the total (pre-biodegradation) electron acceptor concentration 
and F is the mass ratio of EA to ED consumed in the biodegradation reaction.
This equation implies that, concentration being expressed in $mol/L$, that $1mol$ of EA is responsible for degradation of $1mol$ ED.
This is however not the case in most biodegradation reactions, in which case concentrations of ED and EA should be presented in the same stoichiometric units, which is elaborated upon in `Biodegradation stoichiometry`.
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
$C_{ED}(x,y,z,t) = F_{ED}(x,y,z,t, C^0_{ED},\lambda) - (C^0_{EA} - F_{EA}(x,y,z,t,C^0_{EA}))$,
where $F_{ED}$ is a transport equation for a linearly decaying species with source concentration of $C^0_{ED}$ and degradation rate $\lambda$, 
and $F_EA$ is that same transport equation for a conservative species with source concentration $C^0_{EA}$.
All other transport parameters are the same in $F_{ED}$ and $F_{EA}$.
For $F$, Gutierrez-Neri et al. (2009) uses the solution of Domenico and Robbins (1985), with linear degradations as in Domenico (1987).
Hunkeler et al. (2010) pointed out a few issues in the overall formulation of the core-fringe degradation equation, where eq.7 reads:

$$
C_{ED}(x,y,t) = 
\left\{ 
    \begin{array}{l}
        0 \text{ for } C^T_{ED}\cdot \left(K(x,\lambda) \cdot\frac{F_1(x,\lambda,t)}{F_1(x,t)} + \frac{C^0_{EA}}{C^0_{ED}}\right) \leq C_{0,EA}&\\
         C^T_{ED}\cdot \left(K(x,\lambda) \cdot\frac{F_1(x,\lambda,t)}{F_1(x,t)} + \frac{C^0_{EA}}{C^0_{ED}}\right) - C_{0,EA} \text{ elsewhere}
    \end{array}
\right\}
$$
Where

$C^T_{ED} = \frac{C^0_{ED}}{4}\cdot F_1(x,t) \cdot F_2(x,y)$

$K(x,\lambda) = \exp\left[ \left( \frac{x}{2\alpha_x} \right) \left( 1-\sqrt{1+\frac{4\lambda \alpha_x}{v}} \right) \right]$

$F_1(x,t) = \text{erfc}\left(\frac{x-vt}{2\sqrt{\alpha_x-vt}} \right)$

$F_1(x,\lambda,t) = \text{erfc}\left(\frac{x-vt\sqrt{(1+4\alpha_x/v)}}{2\sqrt{\alpha_x-vt}} \right)$

$F_2(x,y) = \left[ \text{erf}\left(\frac{y + Y/2}{2 \sqrt{\alpha_yx}}\right) - \text{erf}\left( \frac{y - Y/2}{2\sqrt{\alpha_yx}} \right) \right]$

