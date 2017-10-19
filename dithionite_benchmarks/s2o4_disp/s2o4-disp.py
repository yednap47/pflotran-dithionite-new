# Id: abiotic.py, Mon 28 Aug 2017 05:03:53 PM MDT #
# Created by Sachin Pandey, Scott Hansen, Satish Karra, LANL
# Description: Benchmark for abiotic reactions.
#------------------------------------------------------------------------------
import sys
import os
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import odespy
sys.path.append('../python')
import pyfun as pf

# ------------------------------------------------------------------------------
# Python numerical simulation
# ------------------------------------------------------------------------------
# Parameters
pars = {'s'				: 1.0, # saturation
		'por'			: 0.25, # porosity
		'v_cell'		: 1.0, # m^3_bulk

		'k_s2o4_disp'	: 1.e-5, # [/s]
}

# Solver options
sopt = {'T' :	10.0 * 24 * 60 * 60, # end of simulation [s from d]
		'N' :	1000, # no of timesteps
}

# Initial conditions
init = {'h'		: 1.0e-11, # [M]
		's2o4'	: 1.e-1, # [M]
		's2o3'	: 1.e-20, # [M]
		'so3'	: 1.e-20, # [M]
		# 'chubbite_vf' : 1-pars['por'], # [m^3/m^3_bulk]
}

dithionite_sandbox = pf.make_dithionite_sandbox(pars)
u, t = pf.run_ode(init, pars, sopt, dithionite_sandbox)

L_water = pars['v_cell'] * pars['por'] * pars['s'] * 1.e3 # [L]
results_ode = {}
results_ode['time'] = t/3600/24 # [d from s]
results_ode['h'] = u[:,0]/L_water
results_ode['s2o4'] = u[:,1]/L_water
results_ode['s2o3'] = u[:,2]/L_water
results_ode['so3'] = u[:,3]/L_water
# results_ode['chubbite_vf'] = u[:,6]

# ------------------------------------------------------------------------------
# Compare with pflotran simulation
# ------------------------------------------------------------------------------
simbasename = "s2o4-disp"
observation_filename = [simbasename + '-obs-0.tec']
variable_list = ['Total H+ [M]', 'Total S2O4-- [M]', 'Total S2O3-- [M]', 'Total SO3-- [M]']
observation_list = ['obs (1)']
results_pflotran =  pf.getobsdata(variable_list=variable_list,observation_list=observation_list,observation_filenames=observation_filename)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
# First plot
fig = plt.figure(figsize=[8,8])
ax = fig.add_subplot(2, 2, 1)
skipfactor = 50 # skip data in ode results
xlims = [0,10]
fontsize = 9
# 
pflo_plotvars = [[variable_list[0]], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars = ['h']
legend_list = ['h - PFLOTRAN','h - odespy']
# 
pf.plot_benchmarks(ax, results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, legend_list=legend_list, xlabel="Time [hr]",ylabel="Concentration [M]", xlims=xlims, skipfactor=skipfactor, fontsize=fontsize)

# Second plot
ax = fig.add_subplot(2, 2, 2)
pflo_plotvars = [[variable_list[1]], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars = ['s2o4']
legend_list = ['s2o4 - PFLOTRAN','s2o4 - odespy']
# 
pf.plot_benchmarks(ax, results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, legend_list=legend_list, xlabel="Time [hr]",ylabel="Concentration [M]", xlims=xlims, skipfactor=skipfactor, fontsize=fontsize)

# Third plot
ax = fig.add_subplot(2, 2, 3)
pflo_plotvars = [[variable_list[2]], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars = ['s2o3']
legend_list = ['s2o3 - PFLOTRAN','s2o3 - odespy']
# 
pf.plot_benchmarks(ax, results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, legend_list=legend_list, xlabel="Time [hr]",ylabel="Concentration [M]", xlims=xlims, skipfactor=skipfactor, fontsize=fontsize)

# Fourth plot
ax = fig.add_subplot(2, 2, 4)
pflo_plotvars = [[variable_list[3]], observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars = ['so3']
legend_list = ['so3 - PFLOTRAN','so3 - odespy']
# 
pf.plot_benchmarks(ax, results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, legend_list=legend_list, xlabel="Time [hr]",ylabel="Concentration [M]", xlims=xlims, skipfactor=skipfactor, fontsize=fontsize)

plt.suptitle("dithionite degradation")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(simbasename + '.png')
