# Transport Model Implementations

Find here the equations that have been implemented in the *mibitrans* package as model classes based on analytical solutions of the ADE under initial and boundary conditions as outlined in the section **Model Setup**.

## *Mibitrans* Model

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

## *Anatrans* Model

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

## *Bioscreen* Model

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

### References

Karanovic, M., C. J. Neville, and C. B. Andrews, *BIOSCREEN-AT: BIOSCREEN with an
exact analytical solution*, Groundwater, 45 (2), 242–245, doi:10.1111/j.1745-6584.2006.00296.x, 2007.

Newell, C. J., R. K. McLeod, and J. R. Gonzales, *BIOSCREEN Natural Attenuation Decision Support System User’s Manual Version 1.3*, Tech. Rep. EPA/600/R-96/087, US EPA, 1996.

Wexler, E. J., *Analytical solutions for one-, two-, and three-dimensional solute transport in ground-water systems with uniform flow*, Tech. Rep. 03-B7, U.S. G.P.O. ; doi:10.3133/twri03B7, 1992
