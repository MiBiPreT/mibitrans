# This example script shows how to make animations of concentration distribution with mibitrans

# The two lines below allow animations to be displayed
import matplotlib

matplotlib.use("TkAgg")

import matplotlib.pyplot as plt
import numpy as np
import mibitrans as mbt

# Define model parameters
hydro = mbt.HydrologicalParameters(
    velocity=0.08,  # Groundwater flow velocity, in [m/day]
    porosity=0.25,  # Effective soil porosity [-]
    alpha_x=5,  # Longitudinal dispersivity, in [m]
    alpha_y=0.2,  # Transverse horizontal dispersivity, in [m]
    alpha_z=0.005,  # Transverse vertical dispersivity, in [m]
)

att = mbt.AttenuationParameters(
    bulk_density=1.7,  # Soil bulk density in [g/m^3]
    partition_coefficient=53,  # Partition coefficient of contaminant to soil organic matter, in [m^3/g]
    fraction_organic_carbon=1e-3,  # Fraction of organic material in the soil [-]
    half_life=0,  # Contaminant half life, in [days]
)

source = mbt.SourceParameters(
    source_zone_boundary=np.array([2, 11, 20]),  # Source zone extent, in [m]
    source_zone_concentration=np.array([15, 4, 1]),  # Source zone concentrations, in [g/m^3]
    depth=3,  # Source depth extent, in [m]
    total_mass=5e6,  # Source mass, in [g] or as np.inf / "infinite"
)

model = mbt.ModelParameters(
    model_length=600,  # Model extent in the longitudinal (x) direction in [m].
    model_width=60,  # Model extent in the transverse horizontal (y) direction in [m].
    model_time=20 * 365,  # Model duration in [days].
    dx=5,  # Model grid discretization step size in the longitudinal (x) direction, in [m].
    dy=0.5,  # Model grid discretization step size in the transverse horizontal (y) direction, in [m].
    dt=182,  # Model time discretization step size, in [days]
)

# Make and run some models, everything shown here also works for anatrans and bioscreen models
mbt_model = mbt.Mibitrans(hydro, att, source, model)
mbt_results_nodecay = mbt_model.run()

mbt_model.attenuation_parameters.half_life = 10 * 365
mbt_results_lindecay = mbt_model.run()

mbt_model.instant_reaction(electron_acceptors=[9, 6, 8, 5, 4], utilization_factor=[3.14, 4.9, 21.8, 4.7, 0.78])
mbt_results_instant = mbt_model.run()


# To produce an animation, set `animate` to True. Also put the animation object to a variable
anim = mbt_results_nodecay.centerline(animate=True)
# Use show to display the animation.
plt.show()

# Combine multiple models by using list input
anim2 = mbt.centerline(
    model=[mbt_results_nodecay, mbt_results_lindecay, mbt_results_instant],
    legend_names=["no decay", "linear decay", "instant reaction"],
    animate=True,
)
plt.show()

# Also works for transverse plots
anim3 = mbt.transverse(
    model=[mbt_results_nodecay, mbt_results_lindecay, mbt_results_instant],
    x_position=150,
    legend_names=["no decay", "linear decay", "instant reaction"],
    animate=True,
)
plt.show()

# For breakthrough curve, animation draws the next point for each time step
anim4 = mbt.breakthrough(
    model=[mbt_results_nodecay, mbt_results_lindecay, mbt_results_instant],
    x_position=150,
    legend_names=["no decay", "linear decay", "instant reaction"],
    animate=True,
)
plt.show()

# Animate 2 and 3d plots as well
anim5 = mbt_results_instant.plume_3d(animate=True, cmap="viridis")
plt.show()

anim6 = mbt_results_nodecay.plume_3d(animate=True, cmap="viridis")
plt.show()

# Animation can be saved to file as .gif, optional argument to adjust amount of frames per second.
anim6.save("plume_animation.gif", fps=10)

# Currently it is not possible to adjust individual line colors/styles when animating multiple models in a single plot.
# A first version of this functionality is found in mibitrans.visualize.animation.animate_1d. However, this is still
# under development and needs improvement to become a core feature of the mibitrans package.
