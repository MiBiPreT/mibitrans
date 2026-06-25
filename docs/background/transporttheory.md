# Transport Theory

## Contaminant Transport Equation

Transport of contaminants in the subsurface is modeled based on the three-dimensional advection–dispersion equation (ADE) with linear equilibrium adsorption for uniform flow in the $x$-direction within a homogeneous isotropic porous medium:

$$
\begin{equation}\tag{1}
\begin{aligned}
R\frac{\partial C}{\partial t}
=
-v\frac{\partial C}{\partial x}
+D_x \frac{\partial^2 C}{\partial x^2}
+D_y \frac{\partial^2 C}{\partial y^2}
+D_z \frac{\partial^2 C}{\partial z^2}
-\mu C
\end{aligned}
\end{equation}
$$

Here $C(x,y,z,t)$ is the contaminant concentration in space $(x \in [0,\infty],y \in [-\infty,\infty],z \in [-\infty,\infty])$ and time $t\geq 0$, $R$ is the linear equilibrium retardation factor, $v$ is the uniform groundwater velocity in $x$-direction. $r_{sinks} = -\mu C$ represents a sink term as result to decay/degradation assuming linear decay at the rate  $\mu$. Note that no decay is naturally included by $\mu = 0$.

$D_{x}$, $D_{y}$ and $D_{y}$ are the longitudinal, transverse horizontal and transverse vertical dispersion coefficients, respectively. Dispersion is modeled following the law of *Scheidegger, [1961]*, simplified to the main flow directions, with $D_i = \alpha_i v$ ($i \in \{x,y,z\}$) where $\alpha_x$, $\alpha_y$ and $\alpha_z$ represent longitudinal, transverse horizontal and transverse vertical dispersivity, respectively.

---

### Boundary Conditions

The contaminant source is modeled as a continuous planar block located at $x=0$ with width $2Y$ and thickness $2Z$:

$$
\begin{equation}\tag{2}
C(0, -Y \lt y \lt Y, -Z \lt z \lt Z; t) = C_s(t)
\end{equation}
$$

Far-field condition:

$$
\begin{equation}\tag{3}
\lim_{\sqrt{x^2+y^2+z^2}\to\infty} C(r,t) = 0
\end{equation}
$$

Zero-gradient conditions at infinity:

$$
\begin{equation}\tag{4}
\frac{\partial C}{\partial x} \to 0, \quad
\frac{\partial C}{\partial y} \to 0, \quad
\frac{\partial C}{\partial z} \to 0
\end{equation}
$$

Initial condition:

$$
\begin{equation} \tag{5}
C(x,y,z,0)=0
\end{equation}
$$

If no source depletion is assumed:

$$
\begin{equation} \tag{6}
C_s(t) = C_0
\end{equation}
$$

---

## Solutions to the Transport Equation

### Exact Solution

*Wexler, [1992]* presented an exact analytical solution to the transport equation (ADE) under the given boundary conditions without adsorption and without source depletion. *West et al, [2007]* reports the solution including linear equilibrium adsorption with:

$$
\begin{equation}\tag{7}
\begin{aligned}
C(x,y,z,t)
&=
\frac{C_0 x}{8\sqrt{\pi \alpha_x v/R}}
\\
&\quad \cdot
\int_0^t
\tau^{-3/2}
\exp\left(
-\mu\tau
-\frac{(x-v\tau/R)^2}{4\alpha_x v\tau/R}
\right)
\\
&\quad \cdot
\left[
\operatorname{erf}
\left(
\frac{y+Y}{2\sqrt{D_y\tau/R}}
\right)
-\operatorname{erf}
\left(
\frac{y-Y}{2\sqrt{D_y\tau/R}}
\right)
\right]
\\
&\quad \cdot
\left[
\operatorname{erf}
\left(
\frac{z+Z}{2\sqrt{D_z\tau/R}}
\right)
-\operatorname{erf}
\left(
\frac{z-Z}{2\sqrt{D_z\tau/R}}
\right)
\right]
d\tau
\end{aligned}
\end{equation}
$$

with $\tau$ being the time integration variable and $\operatorname{erf}$ being the error-function. The solution requires numerical integration which can make its application in screening tools less attractive given the higher computational effort and duration of calculation.

### Integral Approximation

The exact solution of the ADE in Eq.7 can be simplified into a fully analytical form by splitting the effect of time on the different directions following the principle $C(x,y,z,t)/C_0 = C(x,t)\cdot C(y,x)\cdot C(z,x)/C_0$. Specifically, we apply to Eq.7 the substitution $\tau = x/v$ (time traveled of advection front) in the error-function terms for the $y$ and $z$ directions. The remaining integral in $\tau$ in Eq.7 for the transport in $x$-direction can be solved as e.g. shown by *Bear, [1979]*. The resulting approximate transport solution is:

$$
\begin{equation}\tag{8}
\begin{aligned}
C(x,y,z,t)
&=
\frac{C_0}{8}
\left[
\exp\left(\frac{x(1-P)}{2\alpha_x}\right)
\operatorname{erfc}
\left(
\frac{x - P v t/R}
{2\sqrt{\alpha_x v t/R}}
\right)
\right.
\\
&\qquad \quad \ \left. +\exp\left(\frac{x(1+P)}{2\alpha_x}\right)
\operatorname{erfc}
\left(
\frac{x + P v t/R}
{2\sqrt{\alpha_x v t/R}}
\right)
\right]
\\
&\qquad \ \cdot
\left[
\operatorname{erf}
\left(
\frac{y+Y}{2\sqrt{\alpha_y x}}
\right)
-\operatorname{erf}
\left(
\frac{y-Y}{2\sqrt{\alpha_y x}}
\right)
\right]
\\
& \qquad \ \cdot
\left[
\operatorname{erf}
\left(
\frac{z+Z}{2\sqrt{\alpha_z x}}
\right)
-\operatorname{erf}
\left(
\frac{z-Z}{2\sqrt{\alpha_z x}}
\right)
\right]
\end{aligned}
\end{equation}
$$

with

$$
P = \sqrt{1 + \frac{4\mu R\alpha_x}{v}}
$$

The equation provides a closed-form analytical, but approximate solution to the ADE due to the assumptions made for dispersion in the transverse direction. It is thus interesting for screening tool given speed up of calculation, particularly when concentrations in transverse direction are of minor interest.

### Bioscreen Solution

*Bioscreen* *[Newell et al, 1996, 1997]* makes use of the approximate analytical solution of *Domenico, [1987]* in various adaptions. Note that the Domenico-solution has not been derived mathematically rigorous as analytical solution of the ADE, but has been composed of analytical solutions for the individual processes.

The equation of use in *Bioscreen* for the case of linear decay of the contaminant within the plume  is given by:

$$
\begin{equation}\tag{9}
\begin{aligned}
C(x,y,z,t)
&=
\frac{C_0}{8}
\exp
\left[
\frac{x(1-P)}{2\alpha_x}
\right]
\operatorname{erfc}
\left[
\frac{x - P v t/R}
{2\sqrt{\alpha_x v t/R}}
\right]
\\
&\quad \cdot
\left[
\operatorname{erf}
\left(
\frac{y+Y}{2\sqrt{\alpha_y x}}
\right)
-\operatorname{erf}
\left(
\frac{y-Y}{2\sqrt{\alpha_y x}}
\right)
\right]
\\
&\quad \cdot
\left[
\operatorname{erf}
\left(
\frac{z+Z}{2\sqrt{\alpha_z x}}
\right)
-\operatorname{erf}
\left(
\frac{z-Z}{2\sqrt{\alpha_z x}}
\right)
\right]
\end{aligned}
\end{equation}
$$

with

$$
P = \sqrt{1 + \frac{4\mu R\alpha_x}{v}}
$$

The *Bioscreen*-Equation is a simplification of the Integral Solution as it truncates the second term in the sum that describes transport in horizontal direction due to advection and dispersion. Truncating the horizontal transport term was presented by *Sauty, [1980]* in the context of the one-dimensional ADE with constant input. He showed that the approximation is close to the non-truncated solution for Peclet-numbers greater then 100, but deviates for Peclet numbers $\lt 100$ which becomes significant with reducing Peclet numbers. Given that macrodispersion due to aquifer heterogeneity can have significant effect the amount of dispersion, Peclet numbers of $>100$ are not guaranteed in field situation, making the approximation questionable as also discussed by *West et al., [2007]*. Modern options of computation also make it unnecessary.

---

## References

Bear, K., *Hydraulics of groundwater*, London ; New York :: McGraw-Hill International Book Co., 1979.

Domenico, P. A., *An analytical model for multidimensional transport of a decaying contaminant species*, Journal of Hydrology, 91 (1), 49–58, doi:10.1016/0022-1694(87)90127-2, 1987.

Newell, C. J., R. K. McLeod, and J. R. Gonzales, *BIOSCREEN Natural Attenuation Decision Support System User’s Manual Version 1.3*, Tech. Rep. EPA/600/R-96/087, US EPA, 1996.

Newell, C. J., R. K. McLeod, and J. R. Gonzales, *BIOSCREEN Natural Attenuation Decision Support System Version 1.4 Revisions*, Tech. rep., US EPA, 1997.

Sauty, J.-P., *An analysis of hydrodispersive transfer in aquifers*, Water Resour. Res., 16, 145–158, doi:10.1029/WR016i001p00145, 1980.

Scheidegger, A. E., *General theory of dispersion in porous media*, Journal of Geophysical Research, 66(10), 3273-3278, 1961.

West, M. R., B. H. Kueper, and M. J. Ungs, *On the use and error of approximation in
the Domenico (1987) solution*, Groundwater, 45 (2), 126–135, doi:10.1111/j.1745-6584.2006.00280.x, 2007.

Wexler, E. J., *Analytical solutions for one-, two-, and three-dimensional solute transport in ground-water systems with uniform flow*, Tech. Rep. 03-B7, U.S. G.P.O. ; doi:10.3133/twri03B7, 1992
