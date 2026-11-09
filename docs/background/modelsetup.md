# _mibitrans_ Transport Models

## Overview

`mibitrans` implements all three (semi-)analytical solutions to the reactive transport equations (under given boundary conditions) as provided in section *Transport Theory*:

- the `Mibitrans` solutions reflecting the exact solution of *Wexler, 1992* / *West et al, [2007]* 
- the `Anatrans` solution reflecting the *Integral Approximation*
- the `Bioscreen` solution as implemented by *BIOSCREEN* *[Newell et al, 1996]* based on the *Domenico, [1987]* model

All three implemented models package extend these solutions in various aspects:

* Core-fringe-degradation
* Source superposition
* Source depletion

These adaptations can be applied to each of the transport models. The concepts of the first two adaptations are outlined here, while the last one is outlined in a separate section on the *Core-fringe-degradation*. At the end of this section, the specific equations of the transport models as implemented in `mibitrans` are provided. 

---

## Source Superposition
Source geometry of PHC contamination after spill events typically shows a spatially extended source zone where concentrations decrease towards the margin. This effect can be accounted for in analytical modeling making use of the concept of source superposition. 

The source is split into $n$ zones of width $2Y_i$ and source concentration $C_{0,i}$ ($i$ represents the zone index). Specifically, the spatially extended source patch located at $x=0$ and $-Y \lt y \lt Y$ is split into $n$ symmetric source zones in $y$ direction, where each source zone $i$ has a width $2Y_i$ (centered around $y=0$) and a source strength $C_i$. Consequently $\sum_i Y_i = Y$. The values of $C_i$ are assumed to gradually decrease. We deliberately chose the source width and source thickness to be $2Y$ and $2Z$ centered around $y=0$ and $z=0$, respectively. This is different from other source boundary definitions (using source dimensions of $Y/2 \lt y \lt Y/2$) often found in respective literature.

The linear character of the ADE allows to apply the principle of superposition. The concentration distribution of the overall plume $C(x,y,t)$ can be calculated as sum of the individual concentration calculated for each source zone as outlined below. 

---

## Source Depletion

The source strength $C_s(t)$, as defined in the boundary conditions of the ADE (section **Transport Theory**) accounts for a potential change of source concentration over time. Assuming a source depletion rate (also called source decay rate) $\gamma_s$ and a source concentration of $C_0$ at $t=0$, it is defined as

$$
C_s(t) = C_0 \exp{\left(\gamma_s t\right)}    
$$

In literature, this process is also often called source decay - we prefer depletion to avoid mix up of decay processes within the contaminant plume.

The source depletion coefficient $\gamma_s$ is dependent on the amount of contaminant flux out of the source zone. If an initial source mass $m_\mathrm{s,0}$ is defined and the source concentration $C_{0,i}$ for each source zone $i$, i.e. the contaminant concentration released to the groundwater, then the source depletion coefficient is calculated accordingly:

$$
\gamma_s = \frac{Q \cdot \bar{C_{0}}}{m_\mathrm{s,0}}
$$

where $\bar{C_{0}} = \frac{1}{Y}\sum_i Y_i\cdot C_{0,i}$ is the (weighted) average source concentration, $Q = 4\cdot v \cdot \theta_e \cdot Y \cdot Z$ is the volumetric flow rate (other quantities see the parameters section).

If the source mass is considered to be infinite, $\gamma_s$ resolves to zero, and therefore, no source depletion is considered.

---

## Boundary Conditions with Source Adaptations

The implemented transport models in the *mibitrans* package make use of boundary conditions including the source adaptions. Starting point remains the 3D advection–dispersion equation (ADE) with linear equilibrium adsorption for uniform flow in the $x$-direction within a homogeneous isotropic porous medium:

$$
R\frac{\partial C}{\partial t}
=
-v\frac{\partial C}{\partial x}
+D_x \frac{\partial^2 C}{\partial x^2}
+D_y \frac{\partial^2 C}{\partial y^2}
+D_z \frac{\partial^2 C}{\partial z^2}
-\mu C
$$

We have $n$ zones of width $2Y_i$ and source concentration $C_{0,i}$ in each of the zones $i$. The actual boundary condition for each source zone $i$ reads

$$
C_i(x=0, -Y_i^* \lt y \lt Y_i^*, -Z \lt z \lt Z; t) = C_{0,i}^*\,\exp{\left(\gamma_s t\right)}
$$

where $Y_i^*$ is  the net source width and $C_{0,i}^*$ is the net source concentration defined as:

$$
Y_i^* = \sum_j^i Y_j 
$$

and

$$
C_{0,i}^* = C_{0,i} - C_{0,i+1}
$$

where the latter does not apply to the most outer source zone, where $C_{0,n}^* = C_{0,n}$.

The contaminant concentration $C_i$ for each source zone $i$ as result of transport can be calculated with any of the ADE solutions as presented in section **Transport Theory**. The concentration distribution of the overall plume $C(x,y,t)$ is the sum of the individual concentrations calculated for each source zone:

$$
C(x,y,t) = \sum_i^n C_i(x,y,t)
$$

---

## Transport Model Implementations

Find here the equations that have been implemented in the `mibitrans` package as model classes based on analytical solutions of the ADE under adapted initial and boundary conditions. 

### `Mibitrans` Model

The fully exact solution of the ADE accounting for a non-constant, time-dependent source concentrations has been presented in *Karanovic et al., [2007]* as extension of the solution of *Wexler, [1992]*.  

The equation implemented in the *mibitrans* model class, also accounts for the option of source superposition with different source values $C^*_{0,i}$, $i = 1,\ldots,n$. There is no resolution in $z$ direction implemented, thus transport is always evaluated at $z=0$. The implemented equation reads:

The fully exact ADE solution accounting for time-dependent exponential source depletion is:

$$
\begin{aligned}
C(x,y,t)
&=
\sum_{i=1}^{n}
\left\{
C^*_{0,i}
\frac{x}{8\sqrt{\pi \alpha_x v_R}}
\exp(-\gamma_s t)
\right.
\\
&\quad \cdot
\int_{0}^{t}
\frac{1}{\tau^{3/2}}
\exp
\left(
(\gamma_s-\mu)\tau -
\frac{(x-v_R\tau)^2}{4\alpha_x v_R\tau}
\right)
\\
&\quad \cdot
\left[
\operatorname{erfc}
\left(
\frac{y-Y_i}{2\sqrt{\alpha_y v_R\tau}}
\right)-
\operatorname{erfc}
\left(
\frac{y+Y_i}{2\sqrt{\alpha_y v_R\tau}}
\right)
\right]
\\
&\quad \cdot
\left.
\left[
\operatorname{erfc}
\left(
\frac{-Z}{2\sqrt{\alpha_z v_R\tau}}
\right)-
\operatorname{erfc}
\left(
\frac{Z}{2\sqrt{\alpha_z v_R\tau}}
\right)
\right]
d\tau
\right\}
\end{aligned}
$$

with $v_R = \frac{v}{R}$ being the retarded velocity of contaminant transport with $v = \frac{Ki}{\theta_e}$ being the groundwater flow velocity. The meaning of the other parameters are repeated in the table in **Overview Parameters**.

---

### `Anatrans` Model

The equation implemented in the *anatrans* model class is based on the intregral solution (section Transport Theory) that uses a slightly simplified, closed analytical form. As the *mibitrans* model class, it does not resolve the $z$ direction, but transport is always evaluated at $z=0$ and accounts for source depletion and source superposition. The implemented equation reads:

$$
\begin{aligned}
C(x,y,t)
&=
\sum_{i=1}^{n}
\left\{
\frac{C^*_{0,i}}{8}
\exp \left(
-\gamma_s t
\right)
\right.
\\
&\quad \cdot
\left[
\exp \left(
\frac{x\left(1-\tilde P\right)}{2\alpha_x}
\right)
\cdot
\operatorname{erfc} \left(
\frac{x - v_R t \tilde P}{2\sqrt{\alpha_x v_R t }}
\right)
\right.
\\
&\quad \ \,+
\left.
\exp \left(
\frac{x\left(1+\tilde P\right)}{2\alpha_x}
\right)
\cdot
\operatorname{erfc} \left(
\frac{x + v_R t \tilde P}{2\sqrt{\alpha_x v_R t }}
\right)
\right]
\\
&\quad \cdot
\left[
\operatorname{erf} \left(
\frac{y + Y_i}{2\sqrt{\alpha_y x}}
\right)
\operatorname{erf} \left(\
\frac{y - Y_i}{2\sqrt{\alpha_y x}}
\right)
\right]
\\
&\quad \cdot
\left.
\left[
\operatorname{erf} \left(
\frac{Z}{2\sqrt{\alpha_z x}}
\right)
\operatorname{erf} \left(
\frac{-Z}{2\sqrt{\alpha_z x}}
\right)
\right]
\right\}
\end{aligned}
$$

with

$$
\tilde P = \sqrt{1+4 (\mu-\gamma_s) \alpha_x/v_R}
$$

that accounts for the advective components of the source depletion.

The *anatrans* model equation is identical to equation of the Integral Approximation when substituting $P$ with $\tilde P$ and $C_0$ with $C_0 \exp(\gamma_s t)$ and applying the source superposition.

---

### `Bioscreen` Model

The equation implemented in the *bioscreen* model class is based on the *Bioscreen solution* (section Transport Theory) and follows the same choices as used in the Excel based software BIOSCREEN *[Newell et al., 1996]*.

The model implementation also accounts for the option of source depletion and source superposition. Source depletion is handled by accounting for the loss of concentration in reducing the initial mass relative to the time the plume has traveled. Specifically, the use  by replacing the initial mass in the *Bioscreen solution* with $C_s(x,t) = C_0 \exp{\left( -\gamma_s \left( t-\frac{x}{v}\right)\right)}$. Source superposition is reflected by the different source values $C^*_{0,i}$, $i = 1,\ldots,n$. 

$$
\begin{aligned}
C(x,y,t)
&=
\sum_{i=1}^{n}
\left\{
\frac{C^*_{0,i}}{8}
\exp
\left(
-\gamma_s
\left(
t-\frac{x}{v_R}
\right)
\right)
\right.
\\
&\quad \cdot
\exp
\left(
\frac{x(1-P)}{2\alpha_x}
\right)
\operatorname{erfc}
\left(
\frac{x-P v_R t}{2\sqrt{\alpha_x v_R t}}
\right)
\\
&\quad \cdot
\left[
\operatorname{erf}
\left(
\frac{y+Y_i}{2\sqrt{\alpha_y x}}
\right)
-\operatorname{erf}
\left(
\frac{y-Y_i}{2\sqrt{\alpha_y x}}
\right)
\right]
\\
&\quad \cdot
\left.
\left[
\operatorname{erf}
\left(
\frac{Z}{2\sqrt{\alpha_z x}}
\right)
-\operatorname{erf}
\left(
\frac{-Z}{2\sqrt{\alpha_z x}}
\right)
\right]
\right\}
\end{aligned}
$$

with

$$
P = \sqrt{1 + 4\mu \alpha_x / v_R}
$$

Again, the implementation has no dependence on $z$, indicating that it only evaluates the solution at $z=0$ and does not provide a vertically resolved concentration.


---

## References

Domenico, P. A., *An analytical model for multidimensional transport of a decaying contaminant species*, Journal of Hydrology, 91 (1), 49–58, doi:10.1016/0022-1694(87)90127-2, 1987.

Karanovic, M., C. J. Neville, and C. B. Andrews, *BIOSCREEN-AT: BIOSCREEN with an
exact analytical solution*, Groundwater, 45 (2), 242–245, doi:10.1111/j.1745-6584.2006.00296.x, 2007.

Newell, C. J., R. K. McLeod, and J. R. Gonzales, *BIOSCREEN Natural Attenuation Decision Support System User’s Manual Version 1.3*, Tech. Rep. EPA/600/R-96/087, US EPA, 1996.

West, M. R., B. H. Kueper, and M. J. Ungs, *On the use and error of approximation in
the Domenico (1987) solution*, Groundwater, 45 (2), 126–135, doi:10.1111/j.1745-6584.2006.00280.x, 2007.

Wexler, E. J., *Analytical solutions for one-, two-, and three-dimensional solute transport in ground-water systems with uniform flow*, Tech. Rep. 03-B7, U.S. G.P.O. ; doi:10.3133/twri03B7, 1992
