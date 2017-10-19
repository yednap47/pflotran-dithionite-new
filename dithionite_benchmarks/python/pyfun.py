import numpy as np
import itertools as it
import odespy
import matplotlib.pyplot as plt

liter_b_to_m3_b = 1e3

def make_dithionite_sandbox(pars):
	def dithionite_sandbox(u, t):
		# Concentrations
		o2 = u[1]/(pars['por'] * pars['s'] * pars['v_cell'] * 1000.0) # [M]
		s2o4 = u[4]/(pars['por'] * pars['s'] * pars['v_cell'] * 1000.0) # [M]

		# Constants
		L_water = pars['por'] * pars['s'] * pars['v_cell'] * 1000.0 # L_water from m^3_water
		cnv_mobileImmobile = pars['por'] * pars['s'] * 1000.0 # L_h20/m^3_bulk
		
		# DERIVATIVES [mol/s]
		r_s2o4_disp = pars['k_s2o4_disp'] * s2o4 * L_water
		r_s2o4_o2 = pars['k_s2o4_o2'] * s2o4 * o2 * L_water

		return [r_s2o4_disp + 2.0*r_s2o4_o2, # H+'
				0.0 -r_s2o4_o2, # O2(aq)
				0.0, # CrO4--
				0.0, # Cr+++
				-r_s2o4_disp - r_s2o4_o2, # S2O4--
				0.5*r_s2o4_disp, # S2O3--
				r_s2o4_disp+r_s2o4_o2, # SO3--
				r_s2o4_o2, # SO4--
				0.0, # Fe+++
				0.0, # Fe++
				0.0, # fast_Fe++
				0.0, # slow_Fe++
				]
	
	return dithionite_sandbox

def run_ode(init, pars, sopt, function):
	# solver = odespy.RK4(function)
	solver = odespy.CashKarp(function)
	solver.set_initial_condition([
		# MOBILE SPECIES
		init['H+']     * pars['por'] * pars['s'] * pars['v_cell'] * 1000.0, # 0 moles
		init['O2(aq)'] * pars['por'] * pars['s'] * pars['v_cell'] * 1000.0, # 1 moles
		init['CrO4--'] * pars['por'] * pars['s'] * pars['v_cell'] * 1000.0, # 2 moles
		init['Cr+++']  * pars['por'] * pars['s'] * pars['v_cell'] * 1000.0, # 3 moles
		init['S2O4--'] * pars['por'] * pars['s'] * pars['v_cell'] * 1000.0, # 4 moles
		init['S2O3--'] * pars['por'] * pars['s'] * pars['v_cell'] * 1000.0, # 5 moles
		init['SO3--']  * pars['por'] * pars['s'] * pars['v_cell'] * 1000.0, # 6 moles
		init['SO4--']  * pars['por'] * pars['s'] * pars['v_cell'] * 1000.0, # 7 moles
		init['Fe+++']  * pars['por'] * pars['s'] * pars['v_cell'] * 1000.0, # 8 moles
		init['Fe++']   * pars['por'] * pars['s'] * pars['v_cell'] * 1000.0, # 9 moles

		# IMMOBILE SPECIES
		init['fast_Fe++'], # 10 m^3/m^3_bulk
		init['slow_Fe++'] # 11 m^3/m^3_bulk
		])

	time_points = np.linspace(0,sopt['T'],sopt['N']+1)
	u, t = solver.solve(time_points)

	return u, t

def getobsdata(variable_list=[], observation_list=[], observation_filenames=[]):
	"""
	Get observation data from pflotran tec file
	"""
	combined_dict = {}
	for file in observation_filenames:
		variable = []
		f = open(file, 'r')
		title = f.readline()
		title = title.split(',')
		for i in title:
			variable.append(i.strip('"'))
			data = np.genfromtxt(file, skip_header=1)
			data = data.T.tolist()
			var_values_dict = dict(zip(variable, data))
			combined_dict.update(var_values_dict)

		for key in combined_dict.keys():
			if 'Time' in key:
				time = combined_dict[key]

		combined_var_obs_list = [variable_list, observation_list]
		combined_var_obs_list = list(it.product(*combined_var_obs_list))

	combined_dict_trimmed = {}
	combined_dict_trimmed['time']= time
	for item in combined_var_obs_list:
		for key in combined_dict.keys():
			if item[0] in key and item[1] in key:
				var_new = [v for v in combined_dict[key]]
				combined_dict_trimmed[item[0] + " " + item[1]] = var_new

	return combined_dict_trimmed

def plot_benchmarks(ax,results_ode = {}, results_pflotran = {}, ode_plotvars =[], pflo_plotvars = [], legend_list=[], xlabel='', ylabel='', xlims=[], ylims=[], skipfactor=1, fontsize=10, mycmap=plt.cm.jet(np.linspace(0,1,5)), majorFormatter=plt.matplotlib.ticker.FormatStrFormatter("%0.1e")):
	"""
	Plot data to an axis object
	"""
	lns = []
	ctr = 0
	for item in pflo_plotvars:
		ln, = ax.plot(results_pflotran['time'], results_pflotran[item[0] + " " + item[1]], linestyle='-',c=mycmap[ctr])
		lns.append(ln)
		ctr =+ 1

	ctr = 0
	for item in ode_plotvars:
		ln, =  ax.plot(results_ode['time'][::skipfactor],results_ode[item][::skipfactor],ls=' ',marker = 'o',c=mycmap[ctr])
		lns.append(ln)
		ctr =+ 1

	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	if xlims: ax.set_xlim(xlims)
	if ylims: ax.set_ylim(ylims)
	ax.yaxis.set_major_formatter(majorFormatter)
	ax.legend(lns, legend_list, ncol=1, fancybox=True, shadow=False, prop={'size': str(fontsize)}, loc='best')

	return lns
