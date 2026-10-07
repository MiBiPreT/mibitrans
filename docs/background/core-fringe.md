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
Note: For this equation to be correct, $n_1$, $n_2$ and $n_3$ need to be 1 in biodegradation reaction $n_1 ED + n_2 EA \to n_3 \text{Product}$. See `Biodegradation stoichiometry`.
In case of Borden et al. (1986) specifically, Oxygen was used as electron acceptor.
This method of superposition is only valid if the reaction between the EA and ED is rapid compared to the groundwater flow velocity, where EA concentration is the only limiting factor.
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
