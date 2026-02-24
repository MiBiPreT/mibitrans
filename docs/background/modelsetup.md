# Transport Model Setup

## Overview

The implemented models in the *mibitrans* package extend the solutions provided in the section **Transport Theory** in various aspects:

* Source superposition
* Source depletion
* Instant reaction for degradation 

These adaptions can be applied to each of the transport models. The first two are outlined here, while the last one is outlined in a separate section on the *Instant Reaction Model*.

---

## Source Superposition
Source geometry of PHC contamination after spill events typically shows a spatially extended source zone where concentrations decrease towards the margin. This effect can be accounted for in analytical modeling making use of the concept of source superposition. 

The source is split into $n$ zones of width $2Y_i$ and source concentration $C_{0,i}$ ($i$ represents the zone index). Specifically, the spatially extended source patch located at $x=0$ and $Y \lt y \lt Y$ is split into $n$ symmetric source zones in $y$ direction, where each source zone $i$ has a width $2Y_i$ (centered around $y=0$) and a source strength $C_i$. Consequently $\sum_i Y_i = Y$. The values of $C_i$ are assumed to gradually decrease. We deliberately chose the source width and source thickness to be $2Y$ and $2Z$ centered around $y=0$ and $z=0$, respectively. This is different from other source boundary definitions (using source dimensions of $Y/2 \lt y \lt Y/2$) often found in respective literature.

The linear character of the ADE allows to apply the principle of superposition. The concentration distribution of the overall plume $C(x,y,t)$ can be calculated as sum of the individual concentration calculated for each source zone as outlined below. 

---

## Source Depletion

### Concept

The source strength $C_s(t)$, as defined in the boundary conditions of the ADE (section **Transport Theory**) accounts for a potential change of source concentration over time. Assuming a source depletion rate (also called source decay rate) $\gamma_s$ and a source concentration of $C_0$ at $t=0$, it is defined as

$$
C_s(t) = C_0 \exp{\left(\gamma_s t\right)}    
$$

In literature, this process is also often called source decay - we prefer depletion to avoid mix up of decay processes within the contaminant plume.

### Source Depletion Rate

The source depletion coefficient $\gamma_s$ is dependent on the amount of contaminant flux out of the source zone. If an initial source mass $m_\mathrm{s,0}$ is defined and the source concentration $C_{0,i}$ for each source zone $i$, i.e. the contaminant concentration released to the groundwater, then the source depletion coefficient is calculated accordingly:

$$
\gamma_s = \frac{Q \cdot \bar{C_{0}}}{m_\mathrm{s,0}}
$$

where $\bar{C_{0}} = \frac{1}{Y}\sum_i Y_i\cdot C_{0,i}$ is the (weighted) average source concentration, $Q = 4\cdot v \cdot \theta_e \cdot Y \cdot Z$ is the volumetric flow rate (other quantities see the parameters section).

If the source mass is considered to be infinite, $\gamma_s$ resolves to zero, and therefore, no source depletion is considered.

---

## Boundary conditions with Source Adaptions

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
