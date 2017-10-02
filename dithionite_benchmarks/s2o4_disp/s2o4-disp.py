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

		'k_s2o4_disp'	: 0.0, # [/s]
}

# Solver options
sopt = {'T' :	10.0 * 24 * 60 * 60, # end of simulation [s from d]
		'N' :	1000, # no of timesteps
}

# Initial conditions
init = {'C'		: 1.e-2, # [M]
		'D_m'	: 1.e-2, # [M]
		'I'		: 1.e-20, # [M]
		'X'		: 1.e-20, # [M]
		'B'		: 1.e-20, # mol/m^3_bulk
		'D_i'	: 1.e-20, # mol/m^3_bulk
		'chubbite_vf' : 1-pars['por'], # [m^3/m^3_bulk]
}

chrotran_sandbox = pf.make_chrotran_sandbox(pars)
u, t = pf.run_ode(init, pars, sopt, chrotran_sandbox)

# L_water = pars['v_cell'] * pars['por'] * pars['s'] * 1.e3 # [L]
results_ode = {}
results_ode['time'] = t/3600 # [hr from s]
results_ode['C'] = u[:,0]
results_ode['D_m'] = u[:,1]
results_ode['I'] = u[:,2]
results_ode['X'] = u[:,3]
results_ode['B'] = u[:,4]
results_ode['D_i'] = u[:,5]
results_ode['chubbite_vf'] = u[:,6]

# ------------------------------------------------------------------------------
# Compare with pflotran simulation
# ------------------------------------------------------------------------------
simbasename = "abiotic"
observation_filename = [simbasename + '-obs-0.tec']
variable_list = ['Total molasses [M]', 'Total Cr(VI)']
observation_list = ['obs1']
results_pflotran =  pf.getobsdata(variable_list=variable_list,observation_list=observation_list,observation_filenames=observation_filename)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
# First plot
fig = plt.figure(figsize=[5,5])
ax = fig.add_subplot(1, 1, 1)
xlims = [0,0.1]
skipfactor = 50 # skip data in ode results
fontsize = 9

pflo_plotvars = [variable_list, observation_list]
pflo_plotvars = list(it.product(*pflo_plotvars))
ode_plotvars = ['D_m','C']
legend_list = ['D_m - PFLOTRAN','Cr(VI) - PFLOTRAN', 'D_m - odespy','Cr(VI) - odespy']

pf.plot_benchmarks(ax, results_ode=results_ode, results_pflotran=results_pflotran, ode_plotvars=ode_plotvars, pflo_plotvars=pflo_plotvars, legend_list=legend_list, xlabel="Time [hr]",ylabel="Concentration [M]", xlims=xlims, skipfactor=skipfactor, fontsize=fontsize)

plt.suptitle("abiotic benchmark")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(simbasename + '.png')
