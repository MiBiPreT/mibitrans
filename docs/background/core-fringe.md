## Core-fringe degradation model

The concept of a combined core and fringe degradation model is first introduced in Gutierrez-Neri et al. (2009).
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
In case of Borden et al. (1986) specifically, Oxygen was used as electron acceptor and instead of concentrations, they consider ED and EA 'loadings'. 
This method of superposition is only valid if the reaction between the EA and ED is rapid compared to the groundwater flow velocity, where EA concentration is the only limiting factor.
Under this assumption, the biodegradation reaction as an 'instantaneous' reaction, where EA and ED cannot co-exist.
