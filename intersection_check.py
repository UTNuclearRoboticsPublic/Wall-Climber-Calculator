import numpy as np
import math

# L0 : [x,y,z] position of the point
#  v : unit vector pointing towards the intersection point in question
#  r : radius of cylinder centered at [0,0,0] and pointing in x direction
def intersection_check(L0, v, r):
	h = np.array([1,0,0]).reshape(3,1)
	w = L0

	a = v.dot(v) - v.dot(h)**2
	b = 2*(v.dot(w) - (v.dot(h))*(w.dot(h)))
	c = w.dot(w) - w.dot(h)**2 - r**2
    
	check = b**2 - 4*a*c

	if(check < 0):
		# print("There are no intersections")
		return 0
	elif(check == 0):
		# print("There is one intersection")
		return 0
	else:
		# print("There are two intersections")
		t1 = (-b + math.sqrt(check))/(2*a)
		t2 = (-b - math.sqrt(check))/(2*a)
		# print(t1, t2)
		return t1[0]

# points : set of points to intersection check, these are full 4x4 transform matrices
# norm   : Determines whether this should look for the first intersection in the positive or negative direction of the point's normal
# radius_cask : The radius of the cylinder that will be intersected
def intersection_check_list(points, radius_cask, norm):
	if norm > 0:
		norm = 1
	elif norm < 0:
		norm = -1
	distances = [None]*len(points)
	for n in np.arange(len(points)):
		distances[n] = intersection_check(points[n][0:3,3],norm*points[n][0:3,2],radius_cask)

	return distances