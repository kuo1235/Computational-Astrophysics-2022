import ulula.setups.shocktube as setup_shocktube
import ulula.run as ulula_run
import ulula.simulation as ulula_sim
import ulula.plots as ulula_plt
import matplotlib.pyplot as plt
import numpy as np

# setup the Sod shocktube problem in x-direction
setup = setup_shocktube.SetupSodX()

# specify the hydro schemes
hs = ulula_sim.HydroScheme(reconstruction = 'const', limiter = 'minmod', riemann='hll', time_integration='euler', cfl = 0.8)

# run the simulation
sim = ulula_run.run(setup, hydro_scheme=hs, tmax=0.2, nx=100)

# plot the images
q_plot = ['DN','VX','PR']
ulula_plt.plot1d(sim, q_plot)
plt.savefig("shocktube_img.png")
