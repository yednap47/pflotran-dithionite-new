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
# Compare with pflotran simulation
# ------------------------------------------------------------------------------
simbasename = "v2_test2"
observation_filename = [simbasename + '-obs-0.tec']
variable_list = ['Total H+ [M]', 'Total CrO4-- [M]', 'Total O2(aq) [M]', 'Total Fe++ [M]', 'Total Fe+++ [M]', 'Fe(OH)3(s) VF', 'Cr(OH)3(s) VF']
observation_list = ['obs (1)']
results_pflotran =  pf.getobsdata(variable_list=variable_list,observation_list=observation_list,observation_filenames=observation_filename)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
# First plot
fig = plt.figure(figsize=[7,6])
majorFormatter=plt.matplotlib.ticker.FormatStrFormatter("%0.1e")
fontsize = 9

for i in range(0,len(variable_list)):
  ax = fig.add_subplot(3, 3, i+1)
  plt.plot(results_pflotran['time'],results_pflotran[variable_list[i] + ' ' + observation_list[0]])
  plt.title(variable_list[i])
  ax.yaxis.set_major_formatter(majorFormatter)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
plt.savefig(simbasename + '.png')

# ------------------------------------------------------------------------------
# Calculate moles of consumption/productions
# ------------------------------------------------------------------------------
# phi = 0.15
# sat = 1
# results = dict()
# for var in variable_list:
#   if "[mol/m^3]" in var: 
#     results[var.split("[")[0]] = max(results_pflotran[var + ' ' + observation_list[0]])-min(results_pflotran[var + ' ' + observation_list[0]])
#   elif "[M]" in var:
#     results[var.split("[")[0]] = (max(results_pflotran[var + ' ' + observation_list[0]])-min(results_pflotran[var + ' ' + observation_list[0]])) * phi * sat * 1000.0
# 
# # Calculate moles of Fe in Fe(OH)3(s) VF consumed
# rho_rock = 0.0010308000000000123
# mv_fe3 = 34.3600 / 100.0**3 # m^3/mol from cm^3/mol
# fe3_vf = max(results_pflotran['Fe(OH)3(s) VF' + ' ' + observation_list[0]])-min(results_pflotran['Fe(OH)3(s) VF' + ' ' + observation_list[0]])
# results["Fe(OH)3(s)"] = fe3_vf / mv_fe3
# 
# print(results)
