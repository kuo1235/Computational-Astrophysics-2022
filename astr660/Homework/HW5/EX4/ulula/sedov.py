import ulula.examples.examples
import ulula.setups.sedov_taylor as setup_sedov
import ulula.run as ulula_run
import ulula.simulation as ulula_sim
import ulula.plots as ulula_plt
import matplotlib.pyplot as plt
import numpy as np

# setup the Sod shocktube problem in x-direction
setup = setup_sedov.SetupSedov()

# specify the hydro schemes
hs = ulula_sim.HydroScheme(reconstruction = 'const', limiter = 'minmod', riemann='hll', time_integration='euler', cfl = 0.8)

# run the simulation
sim = ulula_run.run(setup, hydro_scheme=hs, tmax=0.02, nx=200)

# plot the 2D images
q_plot = ['DN','PR','VT']
#ulula_plt.plot2d(sim, q_plot=q_plot)
#ulula_plt.plot1d(sim, q_plot=q_plot, plot_type='radius')
ulula_plt.plot1d(sim, q_plot, plot_type='radius',true_solution_func=setup.trueSolution(sim,x,q_plot))
#ulula.examples.examples.sedovTest(nx=200, plot1d=True)
plt.savefig("sedov1_img.png")


# TODO: plot the 2D profiles for density, pressure, and total velocity
