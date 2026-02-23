# Overview Parameters



Parameter used for transport modeling in *mibitrans* package with quantity name, mathematical symbol (used in equations), model input name and default unit. Sorted by parameter type. Quantities marked with a "*" are internally calculated and not defined by the user.

| Parameter                          | Symbol   | *mibitrans* input         | Unit   |
| ---------------------------------- | -------- | ------------------------- |:------:|
| **_Settings_**                     |          |                           |        |
| Longitudinal position              | x        | -                         | m      |
| Transverse horizontal position     | y        | -                         | m      |
| Time                               | t        | -                         | days   |
| Contaminant concentration          | C(x,y,t) | *                         | g/m³   |
| **_Hydrological parameters_**      |          |                           |        |
| Flow velocity                      | v        | velocity                  | m/d    |
| Hydraulic conductivity             | K        | h_conductivity            | m/d    |
| Hydraulic gradient                 | i        | h_gradient                | m/m    |
| Effective porosity                 | θₑ       | porosity                  | -      |
| Longitudinal dispersivity          | αₓ       | alpha_x                   | m      |
| Transverse horizontal dispersivity | αᵧ       | alpha_y                   | m      |
| Transverse vertical dispersivity   | α_z      | alpha_z                   | m      |
| **_Attenuation parameters_**       |          |                           |        |
| Retardation factor                 | R        | retardation               | -      |
| Soil bulk density                  | ρ_b      | bulk_density              | kg/m³  |
| Partition coefficient              | K_oc     | partition_coefficient     | m³/kg  |
| Fraction organic carbon            | f_oc     | fraction_organic_carbon   | -      |
| First order decay coefficient      | μ        | decay_rate                | days⁻¹ |
| Contaminant half-life              | t_1/2    | half_life                 | days   |
| **_Source parameters_**            |          |                           |        |
| Source zone width (from y = 0)     | Y_i      | source_zone_boundary      | m      |
| Source zone concentration          | C₀,i     | source_zone_concentration | g/m³   |
| Net source zone concentration      | C*₀,i    | *                         | g/m³   |
| Total initial source mass          | m_s,0    | total_mass                | g      |
| Source thickness (from z = 0)      | Z        | depth                     | m      |
| **_Model parameters_**             |          |                           |        |
| Model length                       | -        | model_length              | m      |
| Model width                        | -        | model_width               | m      |
| Model time                         | -        | model_time                | days   |
| Length step size                   | -        | dx                        | m      |
| Width step size                    | -        | dy                        | m      |
| Time step size                     | -        | dt                        | days   |
| **_Instant reaction parameters_**  |          |                           |        |
| Delta oxygen concentration         | C_ΔO₂    | delta_oxygen              | g/m³   |
| Delta nitrate concentration        | C_ΔNO₃⁻  | delta_nitrate             | g/m³   |
| Iron (II) concentration            | C_Fe²⁺   | ferrous_iron              | g/m³   |
| Delta sulfate concentration        | C_ΔSO₄²⁻ | delta_sulfate             | g/m³   |
| Methane concentration              | C_CH₄    | methane                   | g/m³   |
| Utilization factor                 | UF       | utilization_factor        | g/g    |
| Biodegradation capacity            | BC       | *                         | g/m³   |
| **_Other_**                        |          |                           |        |
| Source depletion rate              | γ_s      | *                         | days⁻¹ |
| Retarded flow velocity             | v_R      | *                         | m/d    |
