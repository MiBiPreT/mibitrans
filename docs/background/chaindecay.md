# Chain decay

## Background

For certain contaminants, like chlorinated solvents, it is required to track more than a single degradation process. 
Dechlorination of perchloroethene (PCE) to the (relatively) safe ethene involves consecutive degradation to trichloroethene (TCE), dichloroethene (DCE) and vinyl chloride (VC). All four are considered environmental contaminants and known carcinogens.

$PCE \to TCE \to DCE \to VC \to Ethene$

Various tools have been developed to model the fate of these contaminants, such as BIOCHLOR (*Aziz et al., 2000*). BIOCHLOR is similar to BIOSCREEN (*Newell et al., 1996*) in setup and design, but uses a slighly more accurate analytical solution of the contaminant transport equation (the one implemented as `Anatrans` solution in the `mibirem` package). While BIOCHLOR is designed towards modeling chlorinated solvents, the procedure can be applied to any sequential degradation process and is as such easily included to `mibitrans`.

---

## Transport Solution

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

The general ADE for the $n$th compound in the chain ($n\neq1$, i.e. other than the first) reads:

$$
\begin{equation}\tag{2}
\begin{aligned}
R_n\frac{\partial C_n}{\partial t}
=
-v\frac{\partial C_n}{\partial x}
+D_x \frac{\partial^2 C_n}{\partial x^2}
+D_y \frac{\partial^2 C_n}{\partial y^2}
+D_z \frac{\partial^2 C_n}{\partial z^2} \\
+\omega_{n-1} \mu_{n-1} C_{n-1} - \mu_n C_n
\end{aligned}
\end{equation}
$$

The specific form of the decay term follows from the chain of decay. Take the second contaminant in the chain, $X_2$, being the decay product of $X_1$. Beside it's decay term $-\mu_2 C_2$, its ADE contains a source term $\omega_1 \mu_1 C_1$. Here, $\omega$ represents the molar weight ratio between daughter-parent compounds; $\frac{mw_{X_2}}{mw_{X_1}}$. 
Thus, for a chain of degradation reactions of length $n$, the source-sink term reads $\omega_{n-1} \mu_{n-1} C_{n-1} - \mu_n C_n$.

This results in coupled transport equations for all compounds, i.e. $C_2$ is dependent on $C_1$ and therefore cannot be solved as-is.  *Sun et al. (1999)* presented a transformation procedure to express concentrations of subsequent daughter compounds in the concentration of the parent compound. A more detailed explanation and derivation of this procedure is given in *Aziz et al. (2000)*.

For compound $X_2$, the transformed concentration ($A_2$) reads:

$$
\begin{equation}\tag{3}
\begin{aligned}
A_2 = C_2 + C_1\frac{\omega_1\mu_1}{\mu_1-\mu_2}
\end{aligned}
\end{equation}
$$

And for the following daughter compound, it reads:

$$
\begin{equation}\tag{4}
\begin{aligned}
A_3 = C_3 + C_2\frac{\omega_2\mu_2}{\mu_2-\mu_3} + C_1\frac{\omega_1\omega_2\mu_1\mu_2}{(\mu_1-\mu_3)(\mu_2-\mu_3)}
\end{aligned}
\end{equation}
$$

The general transformed form for compound $X_n$ reads:

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

The boundary value for $x=0$ (source zone concentrations) needs similar transformation:

$$
\begin{equation}\tag{7}
\begin{aligned}
A_{0,n} = C_{0,n}\sum^{n-1}_{i=1}C_{0,i}\prod^{n-1}_{m=i}\frac{\omega_m\mu_m}{\mu_m-\mu_n}
\end{aligned}
\end{equation}
$$

Consequently, transport equations of all transformed compounds are uncoupled. They can be solved for $A_n$ as usual, using the analytical solution of choice (`Mibitrans`, `Anatrans`, `Bioscreen`) to the transport equation.

The actual concentration $C_n$ is determined by inverting the transformation in sequence:

$$
\begin{equation}\tag{8}
\begin{aligned}
C_n = A_n - \sum^{n-1}_{i=1}C_i\prod^{n-1}_{m=i}\frac{\omega_m\mu_m}{\mu_m-\mu_n}
\end{aligned}
\end{equation}
$$

---

## Limitations

While retardation can be considered in this transport process, it can only be done so if R is equal for each species. 
Consequently, the transport solutions cannot take difference in transport velocity between species into account. 
A workaround used in BIOCHLOR is taking the median retardation factor of all species.  

Source superposition and source depletion are mutually exclusive for the chain decay solution. 
Currently, `mibitrans` can only handle source superposition for chain decay when source depletion is disabled.

---

## References

Aziz, C. E., Newell, C. J., Gonzales, J. R., Haas, P., Clement, T. P., & Sun, Y. (2000). BIOCHLOR, Natural Attenuation Decision Support System User’s Manual Version 1.0 (p. 54) Manual. EPA.

[Domenico, P., An analytical model for multidimensional transport of a decaying contaminant species, Journal of Hydrology, 91, 49–58, doi:10.1016/0022-1694(87)90127-2, 1987.](https://doi.org/10.1016/0022-1694(87)90127-2)

[Newell, C. J., R. K. Mcleod, J. R. Gonzales, and J. T. Wilson, BIOSCREEN natural attenuation decision support system user’s manual version 1.3, Tech. rep., U.S. EPA, 1996.](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P1007K50.TXT)

[Sun, Y., Petersen, J. N., & Clement, T. P. (1999). Analytical solutions for multiple species reactive transport in multiple dimensions. Journal of Contaminant Hydrology, 35(4), 429–440.](https://doi.org/10.1016/S0169-7722(98)00105-3)
