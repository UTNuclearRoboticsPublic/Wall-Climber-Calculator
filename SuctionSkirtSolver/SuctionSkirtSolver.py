############################################################################################
#      Title     : SuctionSkirtSolver
#      Project   : Wall_Climber_Calculator
#      Copyright : Copyright© The University of Texas at Austin, 2023. All rights reserved.
#                
#          All files within this directory are subject to the following, unless an alternative
#          license is explicitly included within the text of each file.
#
#          This software and documentation constitute an unpublished work
#          and contain valuable trade secrets and proprietary information
#          belonging to the University. None of the foregoing material may be
#          copied or duplicated or disclosed without the express, written
#          permission of the University. THE UNIVERSITY EXPRESSLY DISCLAIMS ANY
#          AND ALL WARRANTIES CONCERNING THIS SOFTWARE AND DOCUMENTATION,
#          INCLUDING ANY WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
#          PARTICULAR PURPOSE, AND WARRANTIES OF PERFORMANCE, AND ANY WARRANTY
#          THAT MIGHT OTHERWISE ARISE FROM COURSE OF DEALING OR USAGE OF TRADE.
#          NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH RESPECT TO THE USE OF
#          THE SOFTWARE OR DOCUMENTATION. Under no circumstances shall the
#          University be liable for incidental, special, indirect, direct or
#          consequential damages or loss of profits, interruption of business,
#          or related expenses which may arise from use of software or documentation,
#          including but not limited to those resulting from defects in software
#          and/or documentation, or loss or inaccuracy of data of any kind.
#
############################################################################################
"""

The code below the functions is designed to calculate through an array of spin
and inversion values for the robot's position and plot them as a surface.

The areas below zero are orientations where the robot will fail by friction
balace, i.e. not be able to drive due to high skirt friction.

The areas highlighted in red are orientations where the robot will fail by
losing suction, i.e. the skirt does not have full contact with
the surface.

"""

import sys
# setting path
sys.path.append('../Wall_Climber_Calculator')
PATH = 'SuctionSkirtSolver/'

import numpy as np
import pandas as pd
from load_yaml import load_yaml
import math
import scipy.optimize as opt
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from friction_balancing_tool import friction_balancing_tool
from plotting_files import main_plotter

# Tests the failure modes at one orientation of the robot
# 	Uses scipy optimize package and the friction balancing tool to descend to the robot position at which
#	forces and moments are balanced. 
#	It then calls the friction balancing tool one more time with outputs set to True and returns the results
#	of failure mode analysis.
def Single_Test(spin, inversion, setup_vars, override_key='', override_var=None, render=False, verbose=False):
	if override_key == "radius_cask":
		radius_cask = override_var
	else:
		radius_cask = setup_vars["radius_cask"]
	
	solution = opt.minimize(friction_balancing_tool, [radius_cask,0], method="Nelder-Mead", args=(spin,inversion,setup_vars,override_key,override_var), bounds=[(radius_cask-5, radius_cask+10),(-0.2, 0.2)])
	results = friction_balancing_tool(solution.x, spin, inversion, setup_vars, override_key, override_var, render=render, verbose=verbose, outputs=True)

	if(verbose):
		print("Roots are:", solution.x)
		print("Success:", solution.success)
		print("Message:", solution.message)
		print("Sliding:", results[2])
		print("Peeling:", results[3])
		if(solution.success):
			print("The Optimization Succeeded to a tolerance of ", results[0], " N and ", results[1], " Nm")
		else:
			print("The Optimization Failed with a tolerance of ", results[0], " N and ", results[1], " Nm")

	return results

# Tests the failure modes across all orientations of the robot
# 	Uses the Single_Test function many times across an array of robot orientations to generate a
#	comprehensive report of each of the failure modes in question. 
#	See plotting_files.py for how to plot and view these results.
def Full_Orientation_Test(resolution,setup_vars,tol=1e-2,override_key='',override_vars=None,plot=False,write=False):
	### Plotting the performance of the robot in an array of orientations ###

	# Only need 0 to pi/2
	spin_max = math.pi/2
	# Only need 0 to pi
	inversion_max = math.pi

	if(write):
		force_columns = [override_key,"spin [rad]","inversion [rad]","force balance"]
		moment_columns = [override_key,"spin [rad]","inversion [rad]","moment balance"]
		friction_columns = [override_key,"spin [rad]","inversion [rad]","friction balance"]
		contact_columns = [override_key,"spin [rad]","inversion [rad]","contact boolean"]
		sliding_columns = [override_key,"spin [rad]","inversion [rad]","sliding balance"]
		force_towrite = []
		moment_towrite = []
		friction_towrite = []
		contact_towrite = []
		sliding_towrite = []

	for override_var in override_vars:
		force_results = np.zeros((resolution,resolution))
		moment_results = np.zeros((resolution,resolution))
		friction_results = np.zeros((resolution,resolution))
		contact_results = np.zeros((resolution,resolution))
		sliding_results = np.zeros((resolution,resolution))
		for i in range(resolution):
			spin = i*spin_max/(resolution-1)
			for j in range(resolution):
				inversion = j*inversion_max/(resolution-1)
				force_results[i,j], moment_results[i,j], friction_results[i,j], contact_results[i,j], sliding_results[i,j] = Single_Test(spin, inversion, setup_vars, override_key=override_key, override_var=override_var)
				force_results[i,j] = 1 if abs(force_results[i,j]) < tol else 0
				moment_results[i,j] = 1 if abs(moment_results[i,j]) < tol else 0
				friction_results[i,j] = friction_results[i,j] if contact_results[i,j] else None
				sliding_results[i,j] = sliding_results[i,j] if contact_results[i,j] else None
				if(write):
					force_towrite.append([override_var, spin, inversion,force_results[i,j]])
					moment_towrite.append([override_var, spin, inversion,moment_results[i,j]])
					friction_towrite.append([override_var, spin, inversion, friction_results[i,j]])
					contact_towrite.append([override_var, spin, inversion, contact_results[i,j]])
					sliding_towrite.append([override_var, spin, inversion, sliding_results[i,j]])
		# if(plot):
		# 	spin = np.linspace(0,resolution,num=resolution,endpoint=True,dtype=int)*spin_max/resolution
		# 	inversion = np.linspace(0,resolution,num=resolution,endpoint=True,dtype=int)*inversion_max/resolution
		# 	plot_tests(spin,inversion,force_results,friction_results,contact_results,sliding_results)

	if(write):
		pd.DataFrame(force_towrite,columns=force_columns).to_csv("results\\force_balance.txt")
		pd.DataFrame(moment_towrite,columns=moment_columns).to_csv("results\\moment_balance.txt")
		pd.DataFrame(friction_towrite,columns=friction_columns).to_csv("results\\friction_balance.txt")
		pd.DataFrame(contact_towrite,columns=contact_columns).to_csv("results\\contact_boolean.txt")
		pd.DataFrame(sliding_towrite,columns=sliding_columns).to_csv("results\\sliding_balance.txt")

# def plot_tests(spin,sliding_results,peeling_results):
# 	resolution_test = len(spin)
# 	tick_res = 5
# 	spin_ticks = list(map(int, np.linspace(0,resolution_test-1,tick_res)))
# 	spin_labels = list(map(round, np.linspace(spin[0],spin[-1],tick_res), map(int,2*np.ones(tick_res))))
# 	plt.subplot(2,1,1)
# 	plt.imshow(sliding_results.T, cmap='RdYlGn_r', vmin = 0, vmax = 1, aspect='auto')
# 	plt.xlabel("Sliding")
# 	plt.xticks(spin_ticks, spin_labels)
# 	plt.subplot(2,1,2)
# 	plt.imshow(peeling_results, cmap='RdYlGn_r', vmin = 0, vmax = 1, aspect='auto')
# 	plt.xlabel("Peeling")
# 	plt.xticks(spin_ticks, spin_labels)
# 	plt.show()

if __name__=='__main__':
	setup_vars = load_yaml(PATH + "setup_variables.yaml")
	radius = setup_vars['radius_cask']
	# resolution = 3
	# Full_Orientation_Test(resolution,setup_vars,override_key='radius_cask',override_vars=radii,write=True)
	# main_plotter()

	spin = math.pi/2
	# Single_Test(spin, setup_vars, override_key='radius_cask', override_var=radii[0], render=True, verbose=True)

	# resolution_test = 100
	# range_test = math.pi/100
	# forcemoment = np.zeros((resolution_test,4))
	# for i in range(resolution_test):
	# 	forcemoment[i,0:3] = forcemoment_balancing_tool([radius-6,-range_test/2+range_test*i/resolution_test], spin, setup_vars, render=False, verbose=False, outputs=True)
	# 	forcemoment[i,3] = forcemoment[i,0] + forcemoment[i,1]
	# 	# forcemoment[i,3] = forcemoment[i,2]**2
	# # print(forcemoment)
	# plt.plot(forcemoment)
	# plt.legend(['Wheel Moment','Array Moment','Grav Moment','Sum of Moments'])
	# plt.show()

	# resolution_test = 1
	# spin_max = math.pi/2
	# sliding = np.zeros(resolution_test)
	# peeling = np.zeros(resolution_test)
	# for i in range(resolution_test):
	# 	temp = Single_Test_Moments(spin_max*i/resolution_test, setup_vars, render=False, verbose=True)
	# 	sliding[i] = temp[2]
	# 	peeling[i] = temp[3]
	# tick_res = 5
	# spin_ticks = list(map(int, np.linspace(0,resolution_test-1,tick_res)))
	# spin_labels = list(map(round, np.linspace(0,spin_max,tick_res), map(int,2*np.ones(tick_res))))
	# plt.subplot(2,1,1)
	# plt.imshow(sliding.reshape((1,resolution_test)), cmap='RdYlGn_r', vmin = 0, vmax = 1, aspect='auto')
	# plt.xlabel("Sliding")
	# plt.xticks(spin_ticks, spin_labels)
	# plt.subplot(2,1,2)
	# plt.imshow(peeling.reshape((1,resolution_test)), cmap='RdYlGn_r', vmin = 0, vmax = 1, aspect='auto')
	# plt.xlabel("Peeling")
	# plt.xticks(spin_ticks, spin_labels)
	# plt.show()

	# print(forcemoment_balancing_tool([2.99310023e+03, 1.92269439e-02], 0, setup_vars, render=True, verbose=True, outputs=True))

	# Single_Test(math.pi/4, 3*math.pi/4, setup_vars, override_key="com_height", override_var=0, render=True, verbose=True)

	override_coms = np.linspace(0,9,10)
	Full_Orientation_Test(3,setup_vars,override_key='com_height',override_vars=override_coms,plot=False,write=True)