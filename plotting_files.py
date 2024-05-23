############################################################################################
#      Title     : plotting_files
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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import pandas as pd
from matplotlib.patches import Rectangle
import mpl_toolkits.mplot3d.art3d as art3d

# Approach to plotting stacked heatmaps borrowed from Paul Brodersen on StackOverflow
def tilted_heatmap_in_3d(arr, z, cmap=cm.RdYlGn, cnorm=colors.TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=1), ax=None):

	alpha_cmap = cmap(np.arange(cmap.N))

	# Set alpha
	alpha_cmap[:,-1] = np.flip(np.linspace(0.25, 0.75, cmap.N, axis=-1))

	# Create new colormap
	alpha_cmap = colors.ListedColormap(alpha_cmap)

	if ax is None:
		fig = plt.figure()
		ax = fig.add_subplot(projection='3d')

	for ii, row in enumerate(arr):
		for jj, value in enumerate(row):
			r = Rectangle((ii-0.5, jj-0.5), 1, 1, color=alpha_cmap(cnorm(value)))
			ax.add_patch(r)
			art3d.pathpatch_2d_to_3d(r, z=z, zdir="z")

	ax.set_xlim(-1, ii+1)
	ax.set_ylim(-1, jj+1)
	ax.set_zlim(0, z)
	ax.get_figure().canvas.draw()

def plot_from_csv(filename, ax, title="No Title"):
	data = pd.read_csv(filename)
	headers = data.columns[1:]
	diameter_set = set(data[headers[0]])

	zmin = min(diameter_set)
	zmax = max(diameter_set)
	spin_max = max(data[headers[1]])
	inversion_max = max(data[headers[2]])
	resolution = 0

	for z in diameter_set:
		subset = data.groupby(data[headers[0]] == z).get_group(True)
		table = subset.pivot_table(index=headers[1],columns=headers[2],values=headers[3])
		arr = np.array(table)
		tilted_heatmap_in_3d(arr,z,ax=ax)
		if(resolution == 0):
			resolution = len(arr)

	tick_res = 3
	spin_ticks = list(map(int, np.linspace(0,resolution-1,tick_res)))
	spin_labels = list(map(round, np.linspace(0,spin_max,tick_res), map(int,2*np.ones(tick_res))))
	inv_ticks = list(map(int, np.linspace(0,resolution-1,tick_res)))
	inv_labels = list(map(round, np.linspace(0,inversion_max,tick_res), map(int,2*np.ones(tick_res))))

	ax.set_xlabel(headers[1])
	ax.set_ylabel(headers[2])
	ax.set_zlabel(headers[0])
	ax.set_zlim(zmin, zmax)
	ax.set_title(title)
	ax.set_xticks(spin_ticks)
	ax.set_xticklabels(spin_labels)
	ax.set_yticks(inv_ticks)
	ax.set_yticklabels(inv_labels)
	return ax

def main_plotter():
	testing_data_forcebalance_filename = "results/Curvature_Testing_Data_forcebalance.csv"
	sim_data_forcebalance_filename = "results/contact_boolean.txt"
	testing_data_frictionbalance_filename = "results/Curvature_Testing_Data_frictionbalance.csv"
	sim_data_frictionbalance_filename = "results/friction_balance.txt"
	sim_success_filename = "results/force_balance.txt"

	fig, ax = plt.subplots(nrows=2,ncols=3,subplot_kw={"projection": "3d"})
	plot_from_csv(sim_success_filename, ax[1,0], "Optimizer Success")
	plot_from_csv(testing_data_forcebalance_filename, ax[0,1], "Full Contact Real")
	plot_from_csv(sim_data_forcebalance_filename, ax[1,1], "Full Contact Sim")
	plot_from_csv(testing_data_frictionbalance_filename, ax[0,2], "Friction Balance Real")
	plot_from_csv(sim_data_frictionbalance_filename, ax[1,2], "Friction Balance Sim")

	plt.show()

if __name__=='__main__':
	main_plotter()