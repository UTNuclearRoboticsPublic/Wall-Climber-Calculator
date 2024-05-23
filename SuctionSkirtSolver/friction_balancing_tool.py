############################################################################################
#      Title     : friction_balancing_tool
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
		phi | rotational position of robot about world x axis  | [rad]
	   spin | rotational position of robot about robot z axis  | [rad]
  inversion | position of robot about cylinder axis			   | [rad]
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
 friction_balance | checks whether the wheel friction is larger than the skirt friction
	 full_contact | checks whether the skirt remains fully sealed
  sliding_balance | checks whether the robot will slide

outputs == False
	  F**2 + M**2 | sum of squares of the forces and moments (to be used by a minimization algorithm)

The first input [h, phi] is designed to be used by an external 
optimizer to find the position of the robot at which all forces
and moments are balanced.

The orientation of the robot is defined by (spin) and (inversion).

(setup_vars) includes many hard coded variables in a dictionary that describe 
the conditions of a rectangular, 4 wheeled, direct drive wall climbing robot.
See setup_variables.py for more details.

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

def friction_balancing_tool(input, spin, inversion, setup_vars, override_key='', override_var=None, render=False, verbose=False, outputs=False):
	[h, phi] = input

	## The easiest way to align this code's spin value with that of my real world data ##
	# This aligns zero spin with the axis of the cask
	spin = spin + math.pi/2

	### Variables ###

	# Orientation
	theta_spin = spin # [rad]
	theta_inversion = inversion # [rad]  #determines the direction of the gravity vector only
	zpos_robot = h # [cm]
	Rspin = rotate(theta_spin, 'z')
	Rinversion = rotate(theta_inversion, 'x')
	Rphi = rotate(phi, 'x')

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
	# Suction
	suction_force = setup_vars['suction_force'] # [N]
	# Skirt
	skirt_resolution = setup_vars['skirt_resolution']
	width_skirt = setup_vars['width_skirt'] # [cm]
	length_skirt = setup_vars['length_skirt'] # [cm]
	thickness_skirt = setup_vars['thickness_skirt'] # [cm]
	zpos_skirt = setup_vars['zpos_skirt'] # [cm]
	stiffness_skirt = setup_vars['stiffness_skirt'] # [N/cm^2]
	static_friction_skirt = setup_vars['static_friction_skirt'] # [~]
	# Wheel
	wheel_resolution = setup_vars['wheel_resolution']
	track_wheel = setup_vars['track_wheel'] # [cm]
	wheelbase_wheel = setup_vars['wheelbase_wheel'] # [cm]
	radius_wheel = setup_vars['radius_wheel'] # [cm]
	height_wheel = setup_vars['height_wheel'] # [cm]
	stiffness_wheel = setup_vars['stiffness_wheel'] # [N/cm]
	static_friction_wheel = setup_vars['static_friction_wheel'] # [~]
	# Cask
	radius_cask = setup_vars['radius_cask'] # [cm]
	height_cask = setup_vars['height_cask'] # [cm]
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

	# Skirt
	skirt_f_T = [None]*skirt_resolution
	skirt_b_T = [None]*skirt_resolution
	skirt_r_T = [None]*skirt_resolution
	skirt_l_T = [None]*skirt_resolution
	for n in np.arange(skirt_resolution):
		skirt_f_T[n] = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(-width_skirt/2+n*width_skirt/(skirt_resolution-1),length_skirt/2,-zpos_skirt)])
		skirt_b_T[n] = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(-width_skirt/2+n*width_skirt/(skirt_resolution-1),-length_skirt/2,-zpos_skirt)])
		skirt_r_T[n] = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(width_skirt/2,-length_skirt/2+(n+1)*length_skirt/(skirt_resolution+1),-zpos_skirt)])
		skirt_l_T[n] = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(-width_skirt/2,-length_skirt/2+(n+1)*length_skirt/(skirt_resolution+1),-zpos_skirt)])
	skirt_T = sum([skirt_f_T, skirt_b_T, skirt_r_T, skirt_l_T],[])
	skirt_C = np.linalg.multi_dot([Rspin,transform(0,0,zpos_robot-zpos_skirt)]) 	

	# External Force
	Fext_T = np.linalg.multi_dot([transform(0,0,zpos_robot),Rphi,Rspin,transform(Fext_c[0],Fext_c[1],Fext_c[2])])

	### Intersection Checks ###
	wheel_t = intersection_check_list(wheel_T, radius_cask, 1)
	skirt_t = intersection_check_list(skirt_T, radius_cask, 1)

	### Calculate Forces ###

	# Gravity
	COM = np.array([0, 0, zpos_robot+com_height])
	G = Rinversion[0:3,0:3].dot(np.array([0, 0, -m*g]))
	if(verbose):
		print("Gravity: ", G[2])

	# Suction
	S = np.array([0, 0, -suction_force])
	if(verbose):
		print("Suction: ", S[2])

	# Wheel Normals & Moments
	wheel_points = np.zeros((4*wheel_resolution,3))
	wheel_F = np.zeros((4*wheel_resolution,3))
	wheel_Normals = np.zeros((4*wheel_resolution))
	wheel_Moments = np.zeros((4*wheel_resolution))
	wheel_Frictions = np.zeros((4*wheel_resolution))
	wheel_out_of_contact = False
	for n in np.arange(len(wheel_T)):
		wheel_Normals[n] = stiffness_wheel*wheel_t[n]*height_wheel/wheel_resolution if wheel_t[n] > 0 else 0 # include wheel force toward surface for convergence
		wheel_points[n] = wheel_T[n][0:3,3].reshape((1,3))
		wheel_r = (wheel_points[n] - wheel_C[0:3,3]).reshape((3,1)) # distance r as a vector
		Normal = wheel_Normals[n] if wheel_t[n] > 0 else 0 # Only provide friction if in contact
		wheel_F[n] = (wheel_T[n][0:3,2]*Normal).reshape((1,3)) # normal force F as a vector
		wheel_Moments[n] = np.dot(np.array([1,0,0]).reshape((1,3)),np.cross(wheel_r,wheel_F[n].reshape((3,1)),axis=0)) # Moment about the y axis, jhat dot (r cross F)
		wheel_Frictions[n] = abs(Normal*static_friction_wheel)

		if wheel_t[n] > 0:
			wheel_out_of_contact = True

	if(verbose):
		print("Wheel Normal: ", sum(wheel_Normals))
		print("Wheel Moment: ", sum(wheel_Moments))
		print("Max Wheel Friction: ", sum(wheel_Frictions))

	# Skirt Normals & Moments
	skirt_points = np.zeros((4*skirt_resolution,3))
	skirt_F = np.zeros((4*skirt_resolution,3))
	skirt_Normals = np.zeros((4*skirt_resolution))
	skirt_Moments = np.zeros((4*skirt_resolution))
	skirt_Frictions = np.zeros((4*skirt_resolution))
	for n in np.arange(len(skirt_T)):
		skirt_Normals[n] = stiffness_skirt*skirt_t[n]*(2*length_skirt + 2*width_skirt)/(4*skirt_resolution) # if skirt_t[n] > 0 else 0
		skirt_points[n] = skirt_T[n][0:3,3].reshape((1,3))
		skirt_r = (skirt_points[n] - skirt_C[0:3,3]).reshape((3,1))
		skirt_F[n] = (skirt_T[n][0:3,2]*skirt_Normals[n]).reshape((1,3))
		skirt_Moments[n] = np.dot(np.array([1,0,0]).reshape((1,3)),np.cross(skirt_r,skirt_F[n].reshape((3,1)),axis=0))		
		skirt_Frictions[n] = abs(skirt_Normals[n]*static_friction_skirt)

	if(verbose):
		print("skirt Normal: ", sum(skirt_Normals))
		print("skirt Moment: ", sum(skirt_Moments))
		print("Max skirt Friction: ", sum(skirt_Frictions))

	# Gravity Moment
	grav_r = (COM.T - wheel_C[0:3,3]).reshape((3,1))
	grav_Moment = np.dot(np.array([1,0,0]).reshape((1,3)),np.cross(grav_r,G,axis=0))[0][0]

	if(verbose):
		print("Gravity Force: ", np.linalg.norm(G))
		print("Gravity Moment: ", grav_Moment)

	# External Force & Moment
	Fext_F = Fext_mag*(Fext_v/np.linalg.norm(Fext_v)) # Normalize unit vector and multiply by magnitude
	Fext_r = (Fext_T[0:3,3] - wheel_C[0:3,3]).reshape((3,1))
	Fext_Moment = np.dot(np.array([0,1,0]).reshape((1,3)),np.cross(Fext_r,Fext_F,axis=0))[0][0]

	if(verbose):
		print("External Force: ", Fext_mag)
		print("External Moment: ", Fext_Moment)

	### Force Balance ###
	F = G[2] + S[2] + sum(wheel_Normals) + sum(skirt_Normals) + Fext_F[2]

	### Moment Balance
	M = sum(wheel_Moments) + sum(skirt_Moments) + grav_Moment + Fext_Moment

	### Calculate Failure Modes ###
	if(outputs):
		# Check Full Contact
		full_contact = True
		# if(any(n<=0 for n in skirt_f_t) or any(n<=0 for n in skirt_b_t) or any(n<=0 for n in skirt_r_t) or any(n<=0 for n in skirt_l_t)):
		if(any(n<=0 for n in skirt_t)):
			full_contact = False

		if(verbose):
			if(full_contact):
				print("The Skirt Has Full Contact")
			else:
				print("The Skirt Loses Full Contact")

		# Friction
		Friction_Balance = sum(wheel_Frictions) - sum(skirt_Frictions)
		friction_balance = 0 if Friction_Balance < 0 else 1
		if(verbose):
			print("Friction Balance: ", Friction_Balance)
			if(friction_balance):
				print("The Robot CAN Move!")
			else:
				print("The Robot CANNOT Move!")

		# Sliding
		Sliding_Balance = sum(wheel_Frictions) + sum(skirt_Frictions) - G[1]
		sliding_balance = 0 if Sliding_Balance < 0 else 1
		if(verbose):
			print("Sliding Balance: ", Sliding_Balance)
			if(sliding_balance):
				print("The robot does NOT slide!")
			else:
				print("The robot does slides!")

	if(render):
		### Rendering Meshes ###
		arrow_scale = 1

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

		# Skirt
		skirt_f = pv.Box(bounds=[skirt_C[0,3]-width_skirt/2, skirt_C[0,3]+width_skirt/2, 
			skirt_C[1,3]+length_skirt/2-thickness_skirt, skirt_C[1,3]+length_skirt/2, 
			skirt_C[2,3]+zpos_skirt, skirt_C[2,3]], level=0, quads=False).rotate_z(theta_spin*180/math.pi, (0,0,0)).rotate_x(phi*180/math.pi, skirt_C[0:3,3])
		skirt_b = pv.Box(bounds=[skirt_C[0,3]-width_skirt/2, skirt_C[0,3]+width_skirt/2, 
			skirt_C[1,3]-length_skirt/2, skirt_C[1,3]-length_skirt/2+thickness_skirt, 
			skirt_C[2,3]+zpos_skirt, skirt_C[2,3]], level=0, quads=False).rotate_z(theta_spin*180/math.pi, (0,0,0)).rotate_x(phi*180/math.pi, skirt_C[0:3,3])
		skirt_r = pv.Box(bounds=[skirt_C[0,3]+width_skirt/2-thickness_skirt, skirt_C[0,3]+width_skirt/2, 
			skirt_C[1,3]-length_skirt/2, skirt_C[1,3]+length_skirt/2, 
			skirt_C[2,3]+zpos_skirt, skirt_C[2,3]], level=0, quads=False).rotate_z(theta_spin*180/math.pi, (0,0,0)).rotate_x(phi*180/math.pi, skirt_C[0:3,3])
		skirt_l = pv.Box(bounds=[skirt_C[0,3]-width_skirt/2, skirt_C[0,3]-width_skirt/2+thickness_skirt, 
			skirt_C[1,3]-length_skirt/2, skirt_C[1,3]+length_skirt/2, 
			skirt_C[2,3]+zpos_skirt, skirt_C[2,3]], level=0, quads=False).rotate_z(theta_spin*180/math.pi, (0,0,0)).rotate_x(phi*180/math.pi, skirt_C[0:3,3])

		skirt_arrows = []
		for n, _ in enumerate(skirt_points):
			skirt_start = np.squeeze(skirt_points[n]).tolist()
			skirt_direction = np.squeeze(skirt_F[n]).tolist()
			skirt_arrows.append(pv.Arrow(start=skirt_start, direction=skirt_direction, scale=np.linalg.norm(skirt_F[n])*arrow_scale))

		# Cask
		cask = pv.Cylinder(radius=radius_cask, height=height_cask, direction=(1,0,0), resolution=100, capping=False).triangulate()

		# Intersection
		wheel_fr_X_cask = wheel_fr.boolean_intersection(cask) 
		wheel_fl_X_cask = wheel_fl.boolean_intersection(cask) 
		wheel_rr_X_cask = wheel_rr.boolean_intersection(cask) 
		wheel_rl_X_cask = wheel_rl.boolean_intersection(cask) 

		skirt_f.flip_normals()
		skirt_b.flip_normals()
		skirt_r.flip_normals()
		skirt_l.flip_normals()
		skirt_f_X_cask = skirt_f.boolean_intersection(cask)
		skirt_b_X_cask = skirt_b.boolean_intersection(cask) 
		skirt_r_X_cask = skirt_r.boolean_intersection(cask) 
		skirt_l_X_cask = skirt_l.boolean_intersection(cask) 

		# Gravity
		gravity_start = np.squeeze(COM).tolist()
		gravity_direction = np.squeeze(G).tolist()
		gravity_arrow = pv.Arrow(start=gravity_start, direction=gravity_direction, scale=np.linalg.norm(G)*arrow_scale)

		# Suction
		suction_start = np.squeeze(COM).tolist()
		suction_direction = np.squeeze(S).tolist()
		suction_arrow = pv.Arrow(start=suction_start, direction=suction_direction, scale=np.linalg.norm(S)*arrow_scale*0.5)

		# Plotting
		p = pv.Plotter()
		p.set_viewup([0,0,1])
		p.set_position([0,100,zpos_robot+100])
		p.set_focus(skirt_C[0:3,3])
		p.add_points(skirt_points)
		p.add_points(wheel_points)
		p.add_mesh(skirt_f, style='wireframe')
		p.add_mesh(skirt_b, style='wireframe')
		p.add_mesh(skirt_r, style='wireframe')
		p.add_mesh(skirt_l, style='wireframe')
		p.add_mesh(wheel_fr, style='wireframe')
		p.add_mesh(wheel_fl, style='wireframe')
		p.add_mesh(wheel_rr, style='wireframe')
		p.add_mesh(wheel_rl, style='wireframe')
		p.add_mesh(gravity_arrow, color=(255,0,0))
		p.add_mesh(suction_arrow, color=(255,0,255))
		for arrow in wheel_arrows:
			p.add_mesh(arrow, color=(0,255,0))
		for arrow in skirt_arrows:
			p.add_mesh(arrow, color=(0,0,255))
		p.add_mesh(cask, style='wireframe')
		try: p.add_mesh(wheel_fr_X_cask)
		except: print("Front Right Wheel No Intersection")
		try: p.add_mesh(wheel_fl_X_cask)
		except: print("Front Left Wheel No Intersection")
		try: p.add_mesh(wheel_rr_X_cask)
		except: print("Rear Right Wheel No Intersection")
		try: p.add_mesh(wheel_rl_X_cask)
		except: print("Rear Left Wheel No Intersection")
		try: p.add_mesh(skirt_f_X_cask)
		except: print("Front Skirt No Intersection")
		try: p.add_mesh(skirt_b_X_cask)
		except: print("Rear Skirt No Intersection")
		try: p.add_mesh(skirt_r_X_cask)
		except: print("Right Skirt No Intersection")
		try: p.add_mesh(skirt_l_X_cask)
		except: print("Left Skirt No Intersection")
		p.show_axes()
		p.show()

	if(outputs):
		return F, M, friction_balance, full_contact, sliding_balance

	return F**2 + M**2 # sum of squares to allow opt.minimize to converge