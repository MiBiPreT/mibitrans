# BIOSCREEN

## General
[*BIOSCREEN*](https://www.epa.gov/water-research/bioscreen-natural-attenuation-decision-support-system) has been developed by the U.S. Environmental Protection Agency (EPA) in collaboration with the U.S. Air Force as natural attenuation decision support tool. Its is meant as screening tool to determine if a full-scale evaluation of a contaminated site is needed. Thereby it is not replacing a more involved numerical model, but serves as preliminary step to evaluate the necessity of an involved numerical model. 

*BIOSCREEN* calculates contaminant concentration distributions in 2D for a constant source under uniform flow conditions based on the advection dispersion equation. The source domain can contain several zones of different concentrations. Three modes of decay can be chose to represent degradation processes:
* no decay
* linear decay
* instantaneous biodegradation reaction

Input parameters for *BIOSCREEN* are few compared to numerical models. After parameter entry, visualization of concentrations at the plume centreline or as a 3D plume can be performed for each of the three decay modes. Version 1.4 provides the option to calculate mass balances of the plume and the source. *BIOSCREEN* comes with two sets of example data.

*BIOSCREEN* is implemented in Excel, and provides a graphical interface. The latest version is *BIOSCREEN* 1.4. It was released to be compatible with Microsoft Excel 5.0. Analysis is performed using macro scripts and cellular calculations. Calculations inside the Excel sheets are hidden in the background, and give insight into the actual calculations behind the model.
The spatial resolution of the plume is fixed to 11 steps over the plume length and 5 steps lateral to the plume. Temporal resolution is fixed to 10 time steps. This resolution is relative to the set model extent and time.

## Differences to `mibitrans`

At short distances, around $x=0$, or early times, *BIOSCREEN* ignores the analytical solution and only uses the source decay term, $\sum_{i=0}^{n} \left[  C^*_{0,i} \exp \left(-k_s \left(t - \frac{xR}{v} \right)\right) \right]$ to calculate concentrations. 

For the *no decay* model, *BIOSCREEN* applies the equation:
\begin{equation}
    C(x, y, t) &= \sum_{i=0}^{n} \Biggl\{ \left( C^*_{0,i} \exp \left[-k_s \left(t - \frac{xR}{v} \right)\right] \right) \bigr. \\
    &\quad \quad \quad \cdot \left\{ \frac{1}{8} \operatorname{erfc} \left[ \frac{x - \frac{vt}{R}}{2\sqrt{\alpha_x \frac{vt}{R}}} \right] \right\}  \\
    &\quad \quad \quad \cdot \left\{ \operatorname{erf} \left[ \frac{y + Y^*_i}{2\sqrt{\alpha_y x}} \right] - \operatorname{erf} \left[ \frac{y - Y^*_i}{2\sqrt{\alpha_y x)}} \right] \right\} \\
    &\quad \quad \quad \cdot \Biggl. \left\{ \operatorname{erf} \left[ \frac{Z}{2\sqrt{\alpha_z x)}} \right] - \operatorname{erf} \left[ \frac{-Z}{2\sqrt{\alpha_z x}} \right] \right\} \Biggr\}
\end{equation}

For the instant reaction model, the equation in *BIOSCREEN* reads:
\begin{equation}
    C(x, y, t) + BC &= \sum_{i=0}^{n} \Biggl\{ \left( C^*_{0,i} \exp \left[-k_s^{inst} \left(t - \frac{xR}{v} \right)\right] + BC \right) \bigr. \\
    &\quad \quad \quad \cdot \left\{ \frac{1}{8} \operatorname{erfc} \left[ \frac{x - \frac{vt}{R}}{2\sqrt{\alpha_x \frac{vt}{R}}} \right] \right\}  \\
    &\quad \quad \quad \cdot \left\{ \operatorname{erf} \left[ \frac{y + Y^*_i}{2\sqrt{\alpha_y x}} \right] - \operatorname{erf} \left[ \frac{y - Y^*_i}{2\sqrt{\alpha_y x)}} \right] \right\} \\
    &\quad \quad \quad \cdot \Biggl. \left\{ \operatorname{erf} \left[ \frac{Z}{2\sqrt{\alpha_z x)}} \right] - \operatorname{erf} \left[ \frac{-Z}{2\sqrt{\alpha_z x}} \right] \right\} \Biggr\} 
\end{equation}

For linear decay, no such additional term is needed to calculate for small distance and times. 

The additional terms $\frac{1}{2} \exp \left[ \frac{xv}{\alpha_xR}\right] \cdot \operatorname{erfc}\left[ \frac{x + \frac{vt}{R}}{2\sqrt{\alpha_x \frac{vt}{R}}}\right]$ in Equation XX are predominantly relevant for small distances and times. These are ignored in *BIOSCREEN*. Consequenly, *BIOSCREEN* results slightly differ from `mibitrans` results for the no decay and instant reaction models. 

## Bugs in *BIOSCREEN*

During the setup of `mibitrans`, *BIOSCREEN* has been thoroughly revised. It was found to be host to a minor erroneous calculation for the instant reaction model. *BIOSCREEN* uses calculations from the *no decay model* and corrects them for the different source decay coefficient and source zone concentrations. However, in these corrections, the wrong source decay coefficient is used, resulting in an underestimation of modelled biodegradation. The size of the error is determined by choice of parameters relating to source decay and biodegradation capacity. 

