# Overview Parameters



Parameter used for transport modeling in *mibitrans* package with quantity name, mathematical symbol (used in equations), model input name and default unit. Sorted by parameter type. Quantities marked with a "*" are internally calculated and not defined by the user.

| Parameter                          | Symbol                 | *mibitrans* input         |    Unit    |
| ---------------------------------- |------------------------| ------------------------- |:----------:|
| **_Settings_**                     |                        |                           |            |
| Longitudinal position              | $x$                    | -                         |    $m$     |
| Transverse horizontal position     | $y$                    | -                         |    $m$     |
| Time                               | $t$                    | -                         |   $day$    |
| Contaminant concentration          | $C(x,y,t)$             | *                         |  $g/m^3$   |
| **_Hydrological parameters_**      |                        |                           |            |
| Flow velocity                      | $v$                    | velocity                  |  $m/day$   |
| Hydraulic conductivity             | $K$                    | h_conductivity            |  $m/day$   |
| Hydraulic gradient                 | $i$                    | h_gradient                |   $m/m$    |
| Effective porosity                 | $θₑ$                   | porosity                  |     -      |
| Longitudinal dispersivity          | $\alpha_x$             | alpha_x                   |    $m$     |
| Transverse horizontal dispersivity | $\alpha_y$             | alpha_y                   |    $m$     |
| Transverse vertical dispersivity   | $\alpha_z$             | alpha_z                   |    $m$     |
| **_Attenuation parameters_**       |                        |                           |            |
| Retardation factor                 | $R$                    | retardation               |     -      |
| Soil bulk density                  | $\rho_b$               | bulk_density              |  $kg/m^3$  |
| Partition coefficient              | $K_{oc}$               | partition_coefficient     |  $m^3/kg$  |
| Fraction organic carbon            | $f_{oc}$               | fraction_organic_carbon   |     -      |
| First order decay coefficient      | $\mu$                  | decay_rate                | $day^{-1}$ |
| Contaminant half-life              | $t_{\frac{1}{2}}$      | half_life                 |   $day$    |
| **_Source parameters_**            |                        |                           |            |
| Source zone width (from y = 0)     | $Y_i$                  | source_zone_boundary      |    $m$     |
| Source zone concentration          | $C_{0,i}$              | source_zone_concentration |  $g/m^3$   |
| Net source zone concentration      | $C_{0,i}^*$            | *                         |  $g/m^3$   |
| Total initial source mass          | $m_{s,0}$              | total_mass                |    $g$     |
| Source thickness (from z = 0)      | $Z$                    | depth                     |    $m$     |
| **_Model parameters_**             |                        |                           |            |
| Model length                       | -                      | model_length              |    $m$     |
| Model width                        | -                      | model_width               |    $m$     |
| Model time                         | -                      | model_time                |   $day$    |
| Length step size                   | -                      | dx                        |    $m$     |
| Width step size                    | -                      | dy                        |    $m$     |
| Time step size                     | -                      | dt                        |   $day$    |
| **_Instant reaction parameters_**  |                        |                           |            |
| Delta oxygen concentration         | $C_{\Delta O_2}$       | delta_oxygen              |  $g/m^3$   |
| Delta nitrate concentration        | $C_{\Delta NO_3^-}$    | delta_nitrate             |  $g/m^3$   |
| Iron (II) concentration            | $C_{\Delta Fe^{2+}}$   | ferrous_iron              |  $g/m^3$   |
| Delta sulfate concentration        | $C_{\Delta SO_4^{2-}}$ | delta_sulfate             |  $g/m^3$   |
| Methane concentration              | $C_{\Delta CH_4}$      | methane                   |  $g/m^3$   |
| Utilization factor                 | $UF$                   | utilization_factor        |   $g/g$    |
| Biodegradation capacity            | $BC$                   | *                         |  $g/m^3$   |
| **_Other_**                        |                        |                           |            |
| Source depletion rate              | $\gamma_s$             | *                         | $day^{-1}$ |
| Retarded flow velocity             | $v_R$                  | *                         |   $m/d$    |
