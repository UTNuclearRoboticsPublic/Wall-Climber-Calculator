############################################################################################
#      Title     : MagneticCrawlerSolver
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
PATH = 'Magnetic/'

import numpy as np
import pandas as pd
from load_yaml import load_yaml
import math
import scipy.optimize as opt
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from forcemoment_balancing_tool import forcemoment_balancing_tool
from plotting_files import main_plotter


def Single_Test_Moments(spin, setup_vars, override_key='', override_var=None, render=False, verbose=False):
	if override_key == "radius_cask":
		radius_cask = override_var
	else:
		radius_cask = setup_vars["radius_cask"]
	
	solution = opt.minimize(forcemoment_balancing_tool, [radius_cask-7,0], method="Nelder-Mead", args=(spin,setup_vars,override_key,override_var), bounds=[(radius_cask-7-5, radius_cask-7+5),(-0.2, 0.2)])
	results = forcemoment_balancing_tool(solution.x, spin, setup_vars, override_key, override_var, render=render, verbose=verbose, outputs=True)

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

# TODO: Fix this
def Full_Spin_Test_Moments(resolution,setup_vars,override_key='',override_vars=None,plot=False,write=False):
	### Plotting the performance of the robot in an array of orientations ###

	# Only need 0 to pi/2
	spin_max = math.pi/2
	override_len = len(override_vars)

	if(write):
		sliding_columns = [override_key,"spin [rad]","sliding"]
		peeling_columns = [override_key,"spin [rad]","peeling"]
		sliding_towrite = []
		peeling_towrite = []

	sliding_results = np.zeros((override_len,resolution))
	peeling_results = np.zeros((override_len,resolution))
	for i, override_var in enumerate(override_vars):
		for j in range(resolution):
			spin = j*spin_max/(resolution-1)
			_, _, sliding_results[i,j], peeling_results[i,j] = Single_Test_Moments(spin, setup_vars, override_key=override_key, override_var=override_var)
			if(write):
				sliding_towrite.append([override_var, spin, sliding_results[i,j]])
				peeling_towrite.append([override_var, spin, peeling_results[i,j]])
	if(plot):
		spin = np.linspace(0,resolution,num=resolution,endpoint=True,dtype=int)*spin_max/resolution
		plot_tests(spin,override_key,override_vars,sliding_results,peeling_results)

	if(write):
		pd.DataFrame(sliding_towrite,columns=sliding_columns).to_csv("results\\sliding.txt")
		pd.DataFrame(peeling_towrite,columns=peeling_columns).to_csv("results\\peeling.txt")

def plot_tests(spin,key,vars,sliding_results,peeling_results):
	resolution_test = len(spin)
	vars_len = len(vars)
	tick_res = 5
	spin_ticks = list(map(int, np.linspace(0,resolution_test-1,tick_res)))
	spin_labels = list(map(round, np.linspace(spin[0],spin[-1],tick_res), map(int,2*np.ones(tick_res))))
	var_ticks = list(map(int, np.linspace(0,vars_len-1,tick_res)))
	var_labels = list(map(round, np.linspace(vars[0],vars[-1],tick_res), map(int,2*np.ones(tick_res))))
	plt.subplot(2,1,1)
	plt.imshow(sliding_results, cmap='RdYlGn_r', vmin = 0, vmax = 1, aspect='auto')
	plt.xlabel("Spin [rad]")
	plt.xticks(spin_ticks, spin_labels)
	plt.ylabel(key)
	plt.yticks(var_ticks, var_labels)
	plt.title("Sliding")
	plt.subplot(2,1,2)
	plt.imshow(peeling_results, cmap='RdYlGn_r', vmin = 0, vmax = 1, aspect='auto')
	plt.xlabel("Spin [rad]")
	plt.xticks(spin_ticks, spin_labels)
	plt.ylabel(key)
	plt.yticks(var_ticks, var_labels)
	plt.title("Peeling")
	plt.show()

if __name__=='__main__':
	setup_vars = load_yaml(PATH + "wombot_variables.yaml")
	# radius = setup_vars['radius_cask']
	# resolution = 3
	# Full_Orientation_Test(resolution,setup_vars,override_key='radius_cask',override_vars=radii,write=True)
	# main_plotter()

	# spin = math.pi/2
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

	# print(forcemoment_balancing_tool([ 1.39301715e+03, -4.56950938e-04], 0, setup_vars, render=True, verbose=True, outputs=True))

	# Single_Test_Moments(0, setup_vars, override_key="zpos_array", override_var=6, render=True, verbose=True)
	Single_Test_Moments(0, setup_vars, render=True, verbose=True)

	# override_coms = np.linspace(0,9,10)
	# override_zpos_array = np.linspace(5.5,6.5,11)
	# override_zpos_array = np.linspace(1,6,6)
	# Full_Spin_Test_Moments(5,setup_vars,override_key='zpos_array',override_vars=override_zpos_array,plot=True,write=False)