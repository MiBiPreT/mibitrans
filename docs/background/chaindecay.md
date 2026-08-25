# Chain decay

## Background

For certain contaminants, like chlorinated solvents, it is required to track more than a single degradation process. 
Dechlorination of perchloroethene (PCE) to the (relatively) save ethene involves consecutive degradation to trichloroethene (TCE), dichloroethene (DCE) and vinyl chloride (VC). 

$PCE \to TCE \to DCE \to VC \to Ethene$

All four of which are considered environmental contaminants and are known carcinogens, with the distribution of VC being the main interest.
Various tools have been developed to model the fate of these contaminants, among which is BIOCHLOR [Aziz et al., 2000].
In setup and design, it is the similar to BIOSCREEN [Newell et al., 1996], aside from including the additional term in the Domenico (1987) solution. (The same as the `anatrans` solution in this package).
While BIOCHLOR is designed towards modeling chlorinated solvents, the procedure can be applied to any sequential degradation process.

## Solution

Consider the governing three-dimensional ADE for degradation of a single contaminant ($X_1$):

$$
\begin{equation}\tag{1}
\begin{aligned}
R_1\frac{\partial C_1}{\partial t}
=
-v\frac{\partial C_1}{\partial x}
+D_x \frac{\partial^2 C_1}{\partial x^2}
+D_y \frac{\partial^2 C_1}{\partial y^2}
+D_z \frac{\partial^2 C_1}{\partial z^2}
-\mu_1 C_1
\end{aligned}
\end{equation}
$$

Here, the sink term $-\mu C_1$ represents first order degradation of the contaminant. 
For the decay product ($X_2$) of this contaminant, the ADE would contain a source term $\omega_1 \mu_1 C_1$, 
here, $\omega$ represents the molar weight ratio between daughter-parent compounds; $\frac{mw_{X_2}}{mw_{X_1}}$.
If compound $X_2$ would subsequently decay, the sink term is similar to $X_1$; $-\mu_2 C_2$. 
Thus, for a chain of degradation reactions of length $n$. The source-sink term would be $\omega_{n-1} \mu_{n-1} C_{n-1} - \mu_n C_n$, for every compound other than the first.
This results in a general ADE for chain decay:

$$
\begin{equation}\tag{2}
\begin{aligned}
R_n\frac{\partial C_n}{\partial t}
=
-v\frac{\partial C_n}{\partial x}
+D_x \frac{\partial^2 C_n}{\partial x^2}
+D_y \frac{\partial^2 C_n}{\partial y^2}
+D_z \frac{\partial^2 C_n}{\partial z^2}
+\omega_{n-1} \mu_{n-1} C_{n-1} - \mu_n C_n
\end{aligned}
\end{equation}
$$

However, these are coupled equations, i.e. $C_2$ is dependent on $C_1$ and therefore cannot be solved as-is. 
In the paper of Sun et al. (1999), a transformation procedure is used to express concentrations of subsequent daughter compounds in the concentration of the parent compound.
A more detailed explanation and derivation of this procedure is given in Aziz et al. (2000) and the aforementioned paper of Sun et al. (1999).
For $X_2$, the transformed concentration ($A_2$) would be:

$$
\begin{equation}\tag{3}
\begin{aligned}
A_2 = C_2 + C_1\frac{\omega_1\mu_1}{\mu_1-\mu_2}
\end{aligned}
\end{equation}
$$

And for the daughter compound of $X_2$, it would be:

$$
\begin{equation}\tag{4}
\begin{aligned}
A_3 = C_3 + C_2\frac{\omega_2\mu_2}{\mu_2-\mu_3} + C_1\frac{\omega_1\omega_2\mu_1\mu_2}{(\mu_1-\mu_3)(\mu_2-\mu_3)}
\end{aligned}
\end{equation}
$$

It follows that for $X_n$, the transformed equation is:

$$
\begin{equation}\tag{5}
\begin{aligned}
A_n = C_n\sum^{n-1}_{i=1}C_i\prod^{n-1}_{m=i}\frac{\omega_m\mu_m}{\mu_m-\mu_n}
\end{aligned}
\end{equation}
$$

Note that for $X_1$, $A_1 = C_1$. With this transformation, the ADE becomes:

$$
\begin{equation}\tag{6}
\begin{aligned}
R_n\frac{\partial A_n}{\partial t}
=
-v\frac{\partial A_n}{\partial x}
+D_x \frac{\partial^2 A_n}{\partial x^2}
+D_y \frac{\partial^2 A_n}{\partial y^2}
+D_z \frac{\partial^2 A_n}{\partial z^2}
-\mu_n A_n
\end{aligned}
\end{equation}
$$

The boundary value for $x=0$ (source zone concentrations) need to be similarly transformed;

$$
\begin{equation}\tag{7}
\begin{aligned}
A_{0,n} = C_{0,n}\sum^{n-1}_{i=1}C_{0,i}\prod^{n-1}_{m=i}\frac{\omega_m\mu_m}{\mu_m-\mu_n}
\end{aligned}
\end{equation}
$$

Now the equations are uncoupled, they can be solved for $A_n$ as usual, using the analytical solution to equation 1 of choice.
Afterward, the actual concentration $C_n$ is determined by inverting the transformation in sequence:

$$
\begin{equation}\tag{8}
\begin{aligned}
C_n = A_n - \sum^{n-1}_{i=1}C_i\prod^{n-1}_{m=i}\frac{\omega_m\mu_m}{\mu_m-\mu_n}
\end{aligned}
\end{equation}
$$

## Limitations

While retardation can be considered in this transport process, it can only be done so if R is equal for each species
i.e. the equation can not take difference in transport velocity between species into consideration.
A 'solution' for this used in BIOCHLOR is to take the median retardation factor of all species.  

Furthermore, source superposition and source depletion are mutually exclusive for the chain decay solution. 
Currently, in this package, only source superposition without source depletion is available for chain decay.

## References

[Aziz, C. E., Newell, C. J., Gonzales, J. R., Haas, P., Clement, T. P., & Sun, Y. (2000). BIOCHLOR, Natural Attenuation Decision Support System User’s Manual Version 1.0 (p. 54) [Manual]. EPA.]

[Domenico, P., An analytical model for multidimensional transport of a decaying contaminant species, Journal of Hydrology, 91, 49–58, doi:10.1016/0022-1694(87)90127-2, 1987.] (https://doi.org/10.1016/0022-1694(87)90127-2)

[Newell, C. J., R. K. Mcleod, J. R. Gonzales, and J. T. Wilson, BIOSCREEN natural attenuation decision support system user’s manual version 1.3, Tech. rep., U.S. EPA, 1996.] https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P1007K50.TXT

[Sun, Y., Petersen, J. N., & Clement, T. P. (1999). Analytical solutions for multiple species reactive transport in multiple dimensions. Journal of Contaminant Hydrology, 35(4), 429–440.] https://doi.org/10.1016/S0169-7722(98)00105-3