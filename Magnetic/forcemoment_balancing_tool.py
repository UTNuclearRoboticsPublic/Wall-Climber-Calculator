############################################################################################
#      Title     : forcemoment_balancing_tool
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

""" Force Calculation Sheet for Wall Climbing Robot on Cylindrical Surface

This function (friction_balancing_tool) is designed to be used by an external
optimizer to find the point at which all forces and moments are balanced.

Inputs

Required
input:=[h, phi]
		  h | distance of robot center from center of cylinder | [cm]
		phi | rotational position of robot about world y axis  | [rad]
	   spin | rotational position of robot about robot z axis  | [rad]
 setup_vars | dictionary input with all required parameters	   | (setup_variables.yaml)

Optional
override_key | dictionary key of setup_var to override 		   | string
override_var | value to replace setup_var with
	  render | set True to render the scene in pyvista
	 verbose | set True to print out variables throughout the calculation
	 outputs | set True to return failure modes instead of force balance value

Output

outputs == True
	  F | sum of forces in normal direction to the cylinder
	  M | sum of moments about the axis of the cylinder
sliding | checks whether the robot will slide
peeling | checks whether the robot will begin to peel off of the wall

outputs == False
	F**2 + M**2 | sum of squares of the forces and moments (to be used by a minimization algorithm)

The first input [h, phi] is designed to be used by an external 
optimizer to find the position of the robot at which all forces
and moments are balanced.

The orientation of the robot is defined by (spin).

(setup_vars) includes many hard coded variables in a dictionary that describe 
the conditions of a rectangular, 4 wheeled, direct drive wall climbing robot.
See wombot_variables.py for more details.

Activating the (render) input variable will render the scene via pyvista in
an interactive window.

Activating the (verbose) input variable will print several important force 
values throughout the calculation, leaving it off will print nothing.

Activating the (outputs) input variable will return a list of useful outputs,
leaving it off will return a force balance value designed for optimization.
"""

import numpy as np
import pyvista as pv
import math
from intersection_check import intersection_check_list

def transform(x, y, z):
	"""Very short and sweet function for creating a 4x4 SE(3) transformation matrix that does not rotate"""
	T = np.eye(4)
	T[0,3] = x
	T[1,3] = y
	T[2,3] = z
	return T

def rotate(theta, axis='x'):
	if axis == 'x':
		R = np.array([[1, 0, 0, 0],
					[0, math.cos(theta), -math.sin(theta), 0],
					[0, math.sin(theta), math.cos(theta), 0],
					[0, 0, 0, 1]]) # about world x axis
		return R
	elif axis == 'y':
		R = np.array([[math.cos(theta), 0, math.sin(theta), 0],
					[0, 1, 0, 0],
					[-math.sin(theta), 0, math.cos(theta), 0],
					[0, 0, 0, 1]]) # about world x axis
		return R
	elif axis == 'z':
		R = np.array([[math.cos(theta), -math.sin(theta), 0, 0],
					[math.sin(theta), math.cos(theta), 0, 0],
					[0, 0, 1, 0],
					[0, 0, 0, 1]]) # about world x axis
		return R
	return 0

# TODO: - remove inversion
#		- remove all unnecessary setup vars
#		- add new setup vars (COM height, magnet resolutions, dimensions, strength)
#		- set gravity direction so that cask is vertical
#		- set render so that cask is vertical
#		- replace array point set with magnet point set
#		- replace array force calc with magnet force calc
#		- remove full contact failure mode
#		- add peeling failure mode analysis - statics equation
#		- inverse the curvature of the cask
def forcemoment_balancing_tool(input, spin, setup_vars, override_key='', override_var=None, render=False, verbose=False, outputs=False):
	[h, phi] = input

	### Variables ###

	# Orientation
	theta_spin = spin # [rad]
	zpos_robot = -h # [cm]
	Rspin = rotate(theta_spin, 'z')
	Rphi = rotate(phi, 'y')

	### Override Var ###
	# this allows the user to methodically vary one variable from the setup
	if override_key in setup_vars.keys():
		setup_vars[override_key] = override_var

	# Loading in the setup variables from the input. 
	# Doing this step here to save time, as the yaml pipeline was added at the last second to keep better track different setups.
	# Gravity
	g = setup_vars['g'] # [m/s^2]
	m = setup_vars['m'] # [kg]
	# Center of Mass
	com_height = setup_vars['com_height'] # [N]
	# Magnet Array
	width_array = setup_vars['width_array'] # [cm]
	width_array_resolution = setup_vars['width_array_resolution']
	length_array = setup_vars['length_array'] # [cm]
	length_array_resolution = setup_vars['length_array_resolution']
	strength_array = setup_vars['strength_array'] # [N*cm^3] Treating this as a linear hooke's law type force for now, will consult with Zack for better modeling
	strength_array_offset = setup_vars['strength_array_offset'] # [cm]
	strength_array_multiplier = setup_vars['strength_array_multiplier'] # Based on chevron tank experiment, the array with backplate supplied more force than individual magnets
	zpos_array = setup_vars['zpos_array'] # [cm]
	# Wheel
	wheel_resolution = setup_vars['wheel_resolution']
	track_wheel = setup_vars['track_wheel'] # [cm]
	wheelbase_wheel = setup_vars['wheelbase_wheel'] # [cm]
	radius_wheel = setup_vars['radius_wheel'] # [cm]
	height_wheel = setup_vars['height_wheel'] # [cm]
	stiffness_wheel = setup_vars['stiffness_wheel'] # [N/cm]
	static_friction_wheel = setup_vars['static_friction_wheel'] # [~]
	# Cask
	height_cask = setup_vars['height_cask'] # [cm]
	radius_cask = setup_vars['radius_cask'] # [cm]
	# External Force
	Fext_c = setup_vars['Fc'] # [cm cm cm] with respect to robot frame
	Fext_v = setup_vars['Fv'] # [cm cm cm] with respect to world frame; unit vector (will be normalized regardless)
	Fext_mag = setup_vars['Fmag'] # [N]

	### Defining Geometries ###
	# Wheels
	theta_wheel = math.atan(wheelbase_wheel/(2*zpos_robot))
	wheel_fr_T = [None]*wheel_resolution
	wheel_fl_T = [None]*wheel_resolution
	wheel_rr_T = [None]*wheel_resolution
	wheel_rl_T = [None]*wheel_resolution
	for n in np.arange(wheel_resolution):
		wheel_fr_T[n] = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(track_wheel/2-height_wheel/2+n*height_wheel/(wheel_resolution-1),wheelbase_wheel/2,0),rotate(-theta_wheel),transform(0,0,-radius_wheel)])
		wheel_fl_T[n] = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(track_wheel/2-height_wheel/2+n*height_wheel/(wheel_resolution-1),-wheelbase_wheel/2,0),rotate(theta_wheel),transform(0,0,-radius_wheel)])
		wheel_rr_T[n] = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(-track_wheel/2-height_wheel/2+n*height_wheel/(wheel_resolution-1),wheelbase_wheel/2,0),rotate(-theta_wheel),transform(0,0,-radius_wheel)])
		wheel_rl_T[n] = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(-track_wheel/2-height_wheel/2+n*height_wheel/(wheel_resolution-1),-wheelbase_wheel/2,0),rotate(theta_wheel),transform(0,0,-radius_wheel)])
	wheel_T = sum([wheel_fr_T, wheel_fl_T, wheel_rr_T, wheel_rl_T],[])
	wheel_C = np.linalg.multi_dot([Rspin,transform(0,0,zpos_robot-radius_wheel)])

	# Magnet Array
	array_T = [None]*width_array_resolution*length_array_resolution
	for m in np.arange(length_array_resolution):
		for n in np.arange(width_array_resolution):
			array_T[m*width_array_resolution+n] = np.linalg.multi_dot([transform(0,0,zpos_robot),
															    Rphi,Rspin,transform(-width_array/2+n*width_array/(width_array_resolution-1),
																-length_array/2+m*length_array/(length_array_resolution-1),
													   			-zpos_array)])
	array_C = np.linalg.multi_dot([Rspin,transform(0,0,zpos_robot)])

	# External Force
	Fext_T = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(Fext_c[0],Fext_c[1],Fext_c[2])])

	### Intersection Checks ###
	wheel_t = intersection_check_list(wheel_T, radius_cask, -1)
	array_t = intersection_check_list(array_T, radius_cask, -1)

	### Calculate Forces ###

	# Gravity
	COM = np.array([0, 0, zpos_robot+com_height])
	G = np.array([-m*g, 0, 0])
	# Gravity Moment
	grav_r = (COM.T - wheel_C[0:3,3]).reshape((3,1))
	grav_Moment = np.dot(np.array([0,1,0]).reshape((1,3)),np.cross(grav_r,G,axis=0))[0][0]

	if(verbose):
		print("Gravity Force: ", np.linalg.norm(G))
		print("Gravity Moment: ", grav_Moment)

	# Wheel Normals & Moments
	wheel_points = np.zeros((4*wheel_resolution,3))
	wheel_F = np.zeros((4*wheel_resolution,3))
	wheel_Normals = np.zeros((4*wheel_resolution))
	wheel_Moments = np.zeros((4*wheel_resolution))
	wheel_Frictions = np.zeros((4*wheel_resolution))
	wheel_out_of_contact = False
	for n in np.arange(len(wheel_T)):
		wheel_Normals[n] = -stiffness_wheel*wheel_t[n]*height_wheel/wheel_resolution #if wheel_t[n] < 0 else 0 # TODO: include wheel force toward surface for convergence
		wheel_points[n] = wheel_T[n][0:3,3].reshape((1,3))
		wheel_r = (wheel_points[n] - wheel_C[0:3,3]).reshape((3,1)) # distance r as a vector
		Normal = wheel_Normals[n] if wheel_t[n] < 0 else 0 # Only provide friction if in contact
		wheel_F[n] = (wheel_T[n][0:3,2]*wheel_Normals[n]).reshape((1,3)) # normal force F as a vector
		wheel_Moments[n] = np.dot(np.array([0,1,0]).reshape((1,3)),np.cross(wheel_r,wheel_F[n].reshape((3,1)),axis=0)) # Moment about the y axis, jhat dot (r cross F)
		wheel_Frictions[n] = abs(Normal*static_friction_wheel)

		if wheel_t[n] > 0:
			wheel_out_of_contact = True

	if(verbose):
		print("Wheel Normal: ", sum(wheel_Normals))
		print("Wheel Moment: ", sum(wheel_Moments))
		print("Max Wheel Friction: ", sum(wheel_Frictions))

	# Magnet Normals & Moments
	array_points = np.zeros((width_array_resolution*length_array_resolution,3))
	array_F = np.zeros((width_array_resolution*length_array_resolution,3))
	array_Normals = np.zeros((width_array_resolution*length_array_resolution))
	array_Moments = np.zeros((width_array_resolution*length_array_resolution))
	for n in np.arange(len(array_T)):
		array_Normals[n] = -strength_array*strength_array_multiplier/(array_t[n]+strength_array_offset)**3 if array_t[n] > 0 else -strength_array*strength_array_multiplier/(strength_array_offset)**3 # inverse cube law
		array_points[n] = array_T[n][0:3,3].reshape((1,3))
		array_r = (array_points[n] - wheel_C[0:3,3]).reshape((3,1))
		array_F[n] = (array_T[n][0:3,2]*array_Normals[n]).reshape((1,3))
		array_Moments[n] = np.dot(np.array([0,1,0]).reshape((1,3)),np.cross(array_r,array_F[n].reshape((3,1)),axis=0))		

	if(verbose):
		print("Array Normal: ", sum(array_Normals))
		print("Array Moment: ", sum(array_Moments))

	# External Force & Moment
	Fext_F = Fext_mag*(Fext_v/np.linalg.norm(Fext_v)) # Normalize unit vector and multiply by magnitude
	Fext_r = (Fext_T[0:3,3] - wheel_C[0:3,3]).reshape((3,1))
	Fext_Moment = np.dot(np.array([0,1,0]).reshape((1,3)),np.cross(Fext_r,Fext_F,axis=0))[0][0]

	if(verbose):
		print("External Force: ", Fext_mag)
		print("External Moment: ", Fext_Moment)

	### Force Balance ###
	F = sum(wheel_Normals) + sum(array_Normals) + Fext_F[2]

	### Moment Balance
	M = sum(wheel_Moments) + sum(array_Moments) + grav_Moment + Fext_Moment

	if(verbose):
		print("Normal Force Balance: ", F)
		print("Peeling Moment Balance: ", M)

	### Calculate Failure Modes ###
	if(outputs):
		## Sliding Condition 
		# if wheel max static friction < gravity then sliding
		sliding = sum(wheel_Frictions) < np.linalg.norm(G)
				
		## Peeling Condition
		# if any wheels are out of contact
		peeling = wheel_out_of_contact
		
		if(verbose):
			print("Sliding failure condition: ", sliding)
			print("Peeling failure condition: ", peeling)

	if(render):
		### Rendering Meshes ###
		arrow_scale = 0.5

		# Wheels
		Ry = rotate(math.pi/2,'y')
		wheel_fr_C = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(track_wheel/2,wheelbase_wheel/2,0),Ry])
		wheel_fl_C = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(track_wheel/2,-wheelbase_wheel/2,0),Ry])
		wheel_rr_C = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(-track_wheel/2,wheelbase_wheel/2,0),Ry])
		wheel_rl_C = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(-track_wheel/2,-wheelbase_wheel/2,0),Ry])
		wheel_fr = pv.Cylinder(radius=radius_wheel, height=height_wheel, center=wheel_fr_C[0:3,3], direction=wheel_fr_C[0:3,0:3].dot(np.array([0,0,1]))).triangulate()
		wheel_fl = pv.Cylinder(radius=radius_wheel, height=height_wheel, center=wheel_fl_C[0:3,3], direction=wheel_fl_C[0:3,0:3].dot(np.array([0,0,1]))).triangulate()
		wheel_rr = pv.Cylinder(radius=radius_wheel, height=height_wheel, center=wheel_rr_C[0:3,3], direction=wheel_rr_C[0:3,0:3].dot(np.array([0,0,1]))).triangulate()
		wheel_rl = pv.Cylinder(radius=radius_wheel, height=height_wheel, center=wheel_rl_C[0:3,3], direction=wheel_rl_C[0:3,0:3].dot(np.array([0,0,1]))).triangulate()
		
		wheel_arrows = []
		for n, _ in enumerate(wheel_points):
			wheel_start = np.squeeze(wheel_points[n]).tolist()
			wheel_direction = np.squeeze(wheel_F[n]).tolist()
			wheel_arrows.append(pv.Arrow(start=wheel_start, direction=wheel_direction, scale=np.linalg.norm(wheel_F[n])*arrow_scale))

		# Magnet Array
		array = pv.Box(bounds=[array_C[0,3]-width_array/2, array_C[0,3]+width_array/2, 
			array_C[1,3]-length_array/2, array_C[1,3]+length_array/2, 
			array_C[2,3]-zpos_array, array_C[2,3]], level=0, quads=False).rotate_z(theta_spin*180/math.pi, (0,0,0)).rotate_y(phi*180/math.pi, array_C[0:3,3])

		array_arrows = []
		for n, _ in enumerate(array_points):
			array_start = np.squeeze(array_points[n]).tolist()
			array_direction = np.squeeze(array_F[n]).tolist()
			array_arrows.append(pv.Arrow(start=array_start, direction=array_direction, scale=np.linalg.norm(array_F[n])*arrow_scale))

		# Cask
		cask_render = pv.Cylinder(radius=radius_cask, height=height_cask, direction=(1,0,0), resolution=1000, capping=False)
		cask = cask_render.triangulate()

		wheel_fr_X_cask = wheel_fr.boolean_difference(cask) 
		wheel_fl_X_cask = wheel_fl.boolean_difference(cask) 
		wheel_rr_X_cask = wheel_rr.boolean_difference(cask) 
		wheel_rl_X_cask = wheel_rl.boolean_difference(cask) 

		array_X_cask = array.boolean_difference(cask) 

		# # Points
		# wheel_points = np.zeros((4*wheel_resolution,3))
		# for n in np.arange(len(wheel_points)):
		# 	wheel_points[n] = wheel_T[n][0:3,3].reshape((1,3))

		# array_points = np.zeros((width_array_resolution*length_array_resolution,3))
		# for n in np.arange(len(array_points)):
		# 	array_points[n] = array_T[n][0:3,3].reshape((1,3))

		# Gravity
		gravity_start = np.squeeze(COM).tolist()
		gravity_direction = np.squeeze(G).tolist()
		gravity_arrow = pv.Arrow(start=gravity_start, direction=gravity_direction, scale=np.linalg.norm(G)*arrow_scale)

		# External Force
		Fext_start = np.squeeze(Fext_T[0:3,3]).tolist()
		Fext_direction = np.squeeze(Fext_v).tolist()
		Fext_arrow = pv.Arrow(start=Fext_start, direction=Fext_direction, scale=Fext_mag*arrow_scale)

		# Plotting
		p = pv.Plotter()
		p.set_viewup([1,0,0])
		p.set_position([0,100,zpos_robot+100])
		p.set_focus(array_C[0:3,3])
		p.add_points(array_points)
		p.add_points(wheel_points)
		p.add_mesh(gravity_arrow, color=(255,0,0))
		p.add_mesh(Fext_arrow, color=(255,255,0))
		p.add_mesh(array, style='wireframe')
		p.add_mesh(wheel_fr, style='wireframe')
		p.add_mesh(wheel_fl, style='wireframe')
		p.add_mesh(wheel_rr, style='wireframe')
		p.add_mesh(wheel_rl, style='wireframe')
		for arrow in wheel_arrows:
			p.add_mesh(arrow, color=(0,255,0))
		for arrow in array_arrows:
			p.add_mesh(arrow, color=(0,0,255))
		p.add_mesh(cask_render)#, style='wireframe')
		try: p.add_mesh(wheel_fr_X_cask)
		except: print("Front Right Wheel No Intersection")
		try: p.add_mesh(wheel_fl_X_cask)
		except: print("Front Left Wheel No Intersection")
		try: p.add_mesh(wheel_rr_X_cask)
		except: print("Rear Right Wheel No Intersection")
		try: p.add_mesh(wheel_rl_X_cask)
		except: print("Rear Left Wheel No Intersection")
		try: p.add_mesh(array_X_cask)
		except: print("Magnet Array No Intersection")
		p.show_axes()
		p.show()

	if(outputs):
		return F, M, sliding, peeling

	return F**2 + M**2 # sum of squares to allow opt.minimize to converge