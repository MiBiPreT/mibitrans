# Mass Balance

Mass balance is calculated based on determined concentrations based on the model choice. At any point in time $t$, the current mass of pure phase contaminant in the source zone is calculated via

$$
m_{\mathrm{s}}(t) = m_\mathrm{s,0}\cdot e^{-\gamma_s \cdot t}
$$

This Equations uses $\gamma_{s}^{inst}$ for the model option *Instant reaction* (find definition of $\gamma_{s}^{inst}$ in the particular section).

The mass of contaminant that has left the source zone is thus:

$$
\Delta m_{\mathrm{s}}(t) = m_\mathrm{s,0} - m_\mathrm{s}(t)
$$

The plume mass is calculated as

$$
m_{\mathrm{p}}(t) = V_{cell} \cdot \sum_{i,j}C(x_{i},y_{j},t)
$$

where summation is over all numerical grid cells defined in $x$- and $y$-direction. $V_{cell}=dx \cdot  dy \cdot  Z \cdot  \theta_e$ is the effective pore volume of the cell with dimensions $dx$, $dy$ and $Z$ being the model resolutions in the $x$ and $y$ direction, and the source depth, respectively. 

The plume mass is dependent on the model dimensions. If the extent of (significant) plume concentrations is larger than the model dimensions, the calculated plume mass does not reflect the actual plume mass. The code implementation checks on this aspect and returns a warning to the user when domain extension is recommended.

In case of no decay during transport ($\mu=0$) mass is conserved and therefore, all mass leaving the source is transported within the contaminant plume. This mass, $m^\mathrm{no \, decay}_{\mathrm{p}} (t)$, can be calculated direct and can thus serve as checkup. The mass that has been transported out of the model zone and the calculated plume mass should be identical:

$$
m^\mathrm{no \, decay}_{\mathrm{p}} (t)= \Delta m_{s}(t)
$$

In case of linear decay or instant reaction, contaminant mass in the plume is reduced due to degradation. Consequently, the mass that left the source zone $\Delta m_{s}(t)$ is the sum of the degraded mass and the plume mass. 

