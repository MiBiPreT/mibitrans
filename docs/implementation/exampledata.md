# Example Data

`mibitrans` comes with three example data sets. A Table below provides an overview of all flow, transport and reaction parameters of the example data sets. Note that values marked with $^*$ follow from relating other parameters, e.g. flow velocity from hydraulic conductivity, porosity and hydraulic gradient or linear decay rate from half-life time.


The standard example data set is chosen independent of a specific field site, but with values that reflect typical field settings. Furthermore, we show-case two field site data sets of the *Bemidji site* and *Keesler site* reflect settings of specific spill events. Find notebooks for these under **Examples**.



| Quantity, Symbol (Unit)                     | Example data | Bemidji site | Keesler site      |
| ------------------------------------------- | ------------ | ------------ | ----------------- |
| Hydraulic conductivity, K (m/d)             | -            | -            | 9.5               |
| Hydraulic gradient, i (m/m)                 | -            | -            | 0.003             |
| Effective porosity, θ_e (-)                 | 0.3          | 0.38         | 0.3               |
| Groundwater flow velocity, v (m/d)          | 0.1          | 0.06         | 0.095*            |
| Longitudinal dispersivity, α_x (m)          | 3            | 0.15         | 9.9               |
| Transverse horizontal dispersivity, α_y (m) | 0.02         | 0.015        | 1.0               |
| Transverse vertical dispersivity, α_z (m)   | 0.005        | 0.0015       | 0                 |
| Soil bulk density, ρ_b (kg/m³)              | 1700         | -            | 1700              |
| Partition coefficient, K_oc (m³/kg)         | 0.053        | -            | 0.038             |
| Fraction organic carbon, f_oc (-)           | 0.001        | -            | 5.7e-5            |
| Retardation coefficient, R (-)              | 1.3*         | 1            | 1.01*             |
| Half-life time (d)                          | 0.5·365      | ∞*           | 0.15·365          |
| Linear decay rate μ (1/d)                   | 0.0038*      | 0            | 0.0127*           |
| Source zone concentrations, C₀,i (g/m³)     | [10, 4, 0.5] | 6            | [13.7, 2.5, 0.06] |
| Source zone boundaries, Y_i (m)             | [5,10,20]    | 1            | [2.1, 11.3, 19.8] |
| Source thickness, Z (m)                     | 2            | 1            | 3.05              |
| Total initial source mass, m_source (g)     | ∞            | ∞            | 2e6               |
| Source decay rate γ_s (1/d)                 | 0*           | 0*           | 1.67e-3*          |
| Available oxygen, C_ΔO₂ (g/m³)              | 9            | 8            | 1.65              |
| Available nitrate, C_ΔNO₃⁻ (g/m³)           | 6            | 0            | 0.07              |
| Available sulfate, C_ΔSO₄²⁻ (g/m³)          | 5            | 0            | 22.4              |
| Iron (II) concentration, C_Fe²⁺ (g/m³)      | 8            | 0            | 16.6              |
| Methane concentration, C_CH₄ (g/m³)         | 4            | 0            | 6.6               |


