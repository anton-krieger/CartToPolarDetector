# This script was created by Anton Krieger, for details on the various 
# calculations you can just think about it yourself :P or just contact me ;)

import numpy as np
import sys
import os

pix2 = 2.*np.pi
pi = np.pi

def get_overall_pixel_position(c):
	""" returns the pixel position in terms of its location in the 
		new coordinate system, e.g., Top left == 'TL' """
	
	if (c[1] >= 0.5):
		if (c[0] >= 0.5): return "TR"
		elif (c[0] <= -0.5): return "TL"
		else: return "T"
	elif (c[1] <= -0.5):
		if (c[0] >= 0.5): return "BR"
		elif (c[0] <= -0.5): return "BL"
		else: return "B"
	else:
		if (c[0] >= 0.5): return "R"
		elif (c[0] <= -0.5): return "L"
		else: return "M"


def get_phi(p):
	""" given a position, it calculates the angle phi """
	
	phi = np.arctan2(p[0],p[1])
	if phi<0.: phi += pix2
	
	return phi

def get_r(p):
	""" given a position, it calculates its vectors length """
	return np.sqrt(p[0]*p[0]+p[1]*p[1])

def get_R_px_max(c,N_x,N_y):
	""" given a position of the new grid in the previous cartesian grid 
		as well as the extend (pixel) of the region, this function
		determines the maxium radius (pixel) that needs to be considered"""
	
	corx = N_x if c[0]<N_x/2. else 0
	cory = N_y if c[1]<N_y/2. else 0
	
	return np.sqrt((c[0]-corx)**2 + (c[1]-cory)**2)

def get_ring_pos(r0,c,R_corner,do):
	""" is applied for the central pixel (pos='M') and determines
		the position of the circle with respect to the four edges of the pixel
		in the end it tells us (in phi direction) which thing is 
		leading to the limit of the area, circle or pixel
		# circle inside pixel: "c"
		# circle outside pixel: "p"
		# first circle, then pixel: "cp"
		# first pixel, then circle: "pc"
		# circle, pixel, circle: "cpc"
	"""

	ring_pos = np.array(["xxx","xxx","xxx","xxx"])
	idx = 0
	if do["has_a_o"] == 0:
		if ((r0 < R_corner[idx%4]) and (r0 < R_corner[(idx-1)%4])): ring_pos[idx] = "c"
		elif ((r0 >= R_corner[idx%4]) and (r0 >= R_corner[(idx-1)%4])): ring_pos[idx] = "p"
	elif do["has_a_o"] == 1:
		if c[0] >= 0: ring_pos[idx] = "pc"
		else: ring_pos[idx] = "cp"
	elif do["has_a_o"] == 2: 
		ring_pos[idx] = "cpc"
	else:
		print("Error! Unexpected number of intersections:",do["has_a_o"])
	
	idx = 1
	if do["has_b_o"] == 0:
		if ((r0 < R_corner[idx%4]) and (r0 < R_corner[(idx-1)%4])): ring_pos[idx] = "c"
		elif ((r0 >= R_corner[idx%4]) and (r0 >= R_corner[(idx-1)%4])): ring_pos[idx] = "p"
	elif do["has_b_o"] == 1:
		if c[1] <= 0: ring_pos[idx] = "pc"
		else: ring_pos[idx] = "cp"
	elif do["has_b_o"] == 2: ring_pos[idx] = "cpc"
	else:
		print("Error! Unexpected number of intersections:",do["has_b_o"])
	
	idx = 2
	if do["has_c_o"] == 0:
		if ((r0 < R_corner[idx%4]) and (r0 < R_corner[(idx-1)%4])): ring_pos[idx] = "c"
		elif ((r0 >= R_corner[idx%4]) and (r0 >= R_corner[(idx-1)%4])): ring_pos[idx] = "p"
	elif do["has_c_o"] == 1:
		if c[0] <= 0: ring_pos[idx] = "pc"
		else: ring_pos[idx] = "cp"
	elif do["has_c_o"] == 2: ring_pos[idx] = "cpc"
	else:
		print("Error! Unexpected number of intersections:",do["has_c_o"])
	
	idx = 3
	if do["has_d_o"] == 0:
		if ((r0 < R_corner[idx%4]) and (r0 < R_corner[(idx-1)%4])): ring_pos[idx] = "c"
		elif ((r0 >= R_corner[idx%4]) and (r0 >= R_corner[(idx-1)%4])): ring_pos[idx] = "p"
	elif do["has_d_o"] == 1:
		if c[1] >= 0: ring_pos[idx] = "pc"
		else: ring_pos[idx] = "cp"
	elif do["has_d_o"] == 2: ring_pos[idx] = "cpc"
	else:
		print("Error! Unexpected number of intersections:",do["has_d_o"])
	
	for string in ring_pos:
		if string not in ["c","p","cp","pc","cpc"]:
			print("Error! Unexpected string: "+string)
			sys.exit()
	
	return ring_pos
	

def get_contr_bins(c,pos,N_phi,N_r,dphi,dr):
	""" the n_phi_min, n_phi_max, n_r_min and n_r_max;
		these are the furthest left, right, bottom, top new pixels getting contributions """

	# 0. preparation: 
	# get four corners in order: TR, BR, BL, TL
	c1 = [c[0] + 0.5, c[1] + 0.5]
	c2 = [c[0] + 0.5, c[1] - 0.5]
	c3 = [c[0] - 0.5, c[1] - 0.5]
	c4 = [c[0] - 0.5, c[1] + 0.5]
	

	# 1. get phi bins with contribution
	if pos == "M": n_phi_min, n_phi_max = 0, N_phi-1
	elif pos == "T":
		n_phi_min = np.floor(get_phi(c3)/dphi)
		n_phi_max = np.floor(get_phi(c2)/dphi)
	else: 
		phis = (get_phi(c1),get_phi(c2),get_phi(c3),get_phi(c4))
		n_phi_min = np.floor(min(ph for ph in phis if ph>0)/dphi)
		n_phi_max = np.floor(max(ph for ph in phis if ph>0)/dphi)
		if pos in ["TL","L"]: 
			# now we check if by concidence, one of the corners is exactly above our defined center (so in phi=0 or phi=2pi direction)
			if 0 in phis: n_phi_max = N_phi-1
		elif pos in ["TR","R"]:
			# now we check if by concidence, one of the corners is exactly above our defined center (so in phi=0 or phi=2pi direction)
			if 0 in phis: n_phi_min = 0
	
	# 2. get r bins with contribution
	rs = (get_r(c1),get_r(c2),get_r(c3),get_r(c4))
	if pos == "M": 
		n_r_min = 0
		n_r_max = min(np.floor(max(rs)/dr),N_r-1)
	elif pos == "T":
		n_r_min = max(np.floor((c[1] - 0.5)/dr),0)
		n_r_max = min(np.floor(max(rs)/dr),N_r-1)
	elif pos == "R":
		n_r_min = max(np.floor((c[0] - 0.5)/dr),0)
		n_r_max = min(np.floor(max(rs)/dr),N_r-1)
	elif pos == "B":
		n_r_min = max(np.floor((-c[1] - 0.5)/dr),0)
		n_r_max = min(np.floor(max(rs)/dr),N_r-1)
	elif pos == "L":
		n_r_min = max(np.floor((-c[0] - 0.5)/dr),0)
		n_r_max = min(np.floor(max(rs)/dr),N_r-1)
	else:
		n_r_min = np.floor(min(rs)/dr)
		n_r_max = min(np.floor(max(rs)/dr),N_r-1)
	
	if n_r_min>N_r-1: n_r_min=-1
	
	return int(n_phi_min), int(n_phi_max), int(n_r_min), int(n_r_max)

def rotate_pixel(c,pos,n_phi_min,n_phi_max,N_phi):
	""" rotates pixel such that its new position would get the
		label 'M', 'T', or 'TR'. The 'M'-pixel is not rotated """
	
	# this is why we require that N_phi%4==0, so that this rotation is
	# possible without issues
	N_4 = int(N_phi/4)
	
	# here the actual rotation is applied
	if pos == "TR": return c,pos,n_phi_min,n_phi_max
	elif pos == "BR":
		c_new = np.array([-c[1],c[0]])
		n_phi_min_new,n_phi_max_new = n_phi_min-N_4,n_phi_max-N_4
		pos_new = "TR"
	elif pos == "BL":
		c_new = np.array([-c[0],-c[1]])
		n_phi_min_new,n_phi_max_new = n_phi_min-2*N_4,n_phi_max-2*N_4
		pos_new = "TR"
	elif pos == "TL":
		c_new = np.array([c[1],-c[0]])
		n_phi_min_new,n_phi_max_new = n_phi_min-3*N_4,n_phi_max-3*N_4
		pos_new = "TR"
	elif pos == "T": return c,pos,n_phi_min,n_phi_max
	elif pos == "R":
		c_new = np.array([-c[1],c[0]])
		n_phi_min_new,n_phi_max_new = n_phi_min+3*N_4,n_phi_max-N_4
		pos_new = "T"
	elif pos == "B":
		c_new = np.array([-c[0],-c[1]])
		n_phi_min_new,n_phi_max_new = n_phi_min+2*N_4,n_phi_max-2*N_4
		pos_new = "T"
	elif pos == "L":
		c_new = np.array([c[1],-c[0]])
		n_phi_min_new,n_phi_max_new = n_phi_min+N_4,n_phi_max-3*N_4
		pos_new = "T"
	elif pos == "M": return c,pos,n_phi_min,n_phi_max
	else:
		print("Error! Unexpected position encountered: ",pos)
		sys.exit()
	
	return c_new,pos_new,n_phi_min_new,n_phi_max_new


def ok(x):
	""" checks whether an intersection occured """
	return (0<=x<=1)

def get_circle_pixel_intersect(c,r0):
	""" - given a circle of radius r0 and a pixel with center c = [c_x,c_y];
		  it calcualtes intersection points of the circle and the pixels edges
		  as well as their phi-angles 
		- returns dictionary with solutions 
	"""
	
	# create dictionary; the o letter symbolizes a circle
	do = {}
	
	# get squared radius
	r2 = r0*r0
	
	# get a0 values
	det = r2 - (c[1] + 0.5)*(c[1] + 0.5)
	if det >= 0:
		
		# get solutions
		t1, det  = -(c[0] - 0.5), np.sqrt(det)
		t1, t2 = t1 - det, t1 + det
		
		# find true intersections, their phi-values and their occurance rates
		if ok(t1): 
			do["a_o"] = t1
			do["phi_a_o"] = np.arctan2(-det,c[1]+0.5)
			if do["phi_a_o"] < 0: do["phi_a_o"] += pix2
			if ok(t2): 
				do["a_o+"] = t2
				do["phi_a_o+"] = np.arctan2(det,c[1]+0.5)
				if do["phi_a_o+"] < 0: do["phi_a_o+"] += pix2
				if do["phi_a_o"] < np.pi: do["phi_a_o"] += pix2
				do["has_a_o"] = 2
			else: do["has_a_o"] = 1
		elif ok(t2): 
			do["a_o"] = t2
			do["phi_a_o"] = np.arctan2(det,c[1]+0.5)
			if do["phi_a_o"] < 0: do["phi_a_o"] += pix2
			do["has_a_o"] = 1
		else:
			do["has_a_o"] = 0
	else: do["has_a_o"] = 0

	
	# get b0 values
	det = r2 - (c[0] + 0.5)*(c[0] + 0.5)
	if det >= 0:
		
		# get solutions
		t1, det  = (c[1] + 0.5), np.sqrt(det)
		t1, t2 = t1 - det, t1 + det
		
		# find true intersections, their phi-values and their occurance rates
		if ok(t1): 
			do["b_o"] = t1
			do["phi_b_o"] = np.arctan2(c[0]+0.5,det)
			if do["phi_b_o"] < 0: do["phi_b_o"] += pix2
			if ok(t2): 
				do["b_o+"] = t2
				do["phi_b_o+"] = np.arctan2(c[0]+0.5,-det)
				if do["phi_b_o+"] < 0: do["phi_b_o+"] += pix2
				do["has_b_o"] = 2
			else: do["has_b_o"] = 1
		elif ok(t2): 
			do["b_o"] = t2
			do["phi_b_o"] = np.arctan2(c[0]+0.5,-det)
			if do["phi_b_o"] < 0: do["phi_b_o"] += pix2
			do["has_b_o"] = 1
		else:
			do["has_b_o"] = 0
	else: do["has_b_o"] = 0

	
	# get c0 values
	det = r2 - (c[1] - 0.5)*(c[1] - 0.5)
	if det >= 0:
		
		# get solutions
		t1, det  = (c[0] + 0.5), np.sqrt(det)
		t1, t2 = t1 - det, t1 + det
		
		# find true intersections, their phi-values and their occurance rates
		if ok(t1): 
			do["c_o"] = t1
			do["phi_c_o"] = np.arctan2(det,c[1]-0.5)
			if do["phi_c_o"] < 0: do["phi_c_o"] += pix2
			if ok(t2): 
				do["c_o+"] = t2
				do["phi_c_o+"] = np.arctan2(-det,c[1]-0.5)
				if do["phi_c_o+"] < 0: do["phi_c_o+"] += pix2
				do["has_c_o"] = 2
			else: do["has_c_o"] = 1
		elif ok(t2): 
			do["c_o"] = t2
			do["phi_c_o"] = np.arctan2(-det,c[1]-0.5)
			if do["phi_c_o"] < 0: do["phi_c_o"] += pix2
			do["has_c_o"] = 1
		else:
			do["has_c_o"] = 0
	else: do["has_c_o"] = 0

	
	# get d0 values
	det = r2 - (c[0] - 0.5)*(c[0] - 0.5)
	if det >= 0:
		
		# get solutions
		t1, det  = -(c[1] - 0.5), np.sqrt(det)
		t1, t2 = t1 - det, t1 + det
		
		# find true intersections, their phi-values and their occurance rates
		if ok(t1): 
			do["d_o"] = t1
			do["phi_d_o"] = np.arctan2(c[0]-0.5,-det)
			if do["phi_d_o"] < 0: do["phi_d_o"] += pix2
			if ok(t2): 
				do["d_o+"] = t2
				do["phi_d_o+"] = np.arctan2(c[0]-0.5,det)
				if do["phi_d_o+"] < 0: do["phi_d_o+"] += pix2
				do["has_d_o"] = 2
			else: do["has_d_o"] = 1
		elif ok(t2): 
			do["d_o"] = t2
			do["phi_d_o"] = np.arctan2(c[0]-0.5,det)
			if do["phi_d_o"] < 0: do["phi_d_o"] += pix2
			do["has_d_o"] = 1
		else:
			do["has_d_o"] = 0
	else: do["has_d_o"] = 0

	return do


def get_ray_pixel_intersect_M(ce,phi,Phi_corner):
	""" get the intersection of a ray and the rim of a pixel 
	    assuming the origin of the ray (center) is inside the pixel """
	dI = {}
	
	if phi<=Phi_corner[0]:
		# we need to find an a value 
		a = (ce[1]+0.5) * np.tan(phi) - (ce[0]-0.5)
		if -1e-12< a <1.+1e-12: 
			# ~ dI["a_I"] = max(min(a,1.),0.)
			dI["a_I"] = a		# we accept a small exceedance of the limits
		else: 
			print("Error! the intersection is to far from the expected region:",a)
			sys.exit()
	elif phi<=Phi_corner[1]:
		# we need to find an b value 
		b = - (ce[0]+0.5) / np.tan(phi) + (ce[1]+0.5)
		if -1e-12< b <1.+1e-12: 
			# ~ dI["b_I"] = max(min(b,1.),0.)
			dI["b_I"] = b		# we accept a small exceedance of the limits
		else: 
			print("Error! the intersection is to far from the expected region:",a)
			sys.exit()
	elif phi<=Phi_corner[2]:
		# we need to find an c value 
		c = - (ce[1]-0.5) * np.tan(phi) + (ce[0]+0.5)
		if -1e-12< c <1.+1e-12: 
			# ~ dI["c_I"] = max(min(c,1.),0.)
			dI["c_I"] = c		# we accept a small exceedance of the limits
		else: 
			print("Error! the intersection is to far from the expected region:",a)
			sys.exit()
	elif phi<=Phi_corner[3]:
		# we need to find an d value 
		d = (ce[0]-0.5) / np.tan(phi) - (ce[1]-0.5)
		if -1e-12< d <1.+1e-12: 
			# ~ dI["d_I"] = max(min(d,1.),0.)
			dI["d_I"] = d		# we accept a small exceedance of the limits
		else: 
			print("Error! the intersection is to far from the expected region:",a)
			sys.exit()
	else:
		# we need to find an a value 
		a = (ce[1]+0.5) * np.tan(phi) - (ce[0]-0.5)
		if -1e-12< a <1.+1e-12: 
			# ~ dI["a_I"] = max(min(a,1.),0.)
			dI["a_I"] = a		# we accept a small exceedance of the limits
		else: 
			print("Error! the intersection is to far from the expected region:",a)
			sys.exit()
	
	return dI
		
def get_ray_pixel_intersect_non_M(ce,phi):
	""" - this routine uses the phi-value of a ray and the position
	   information of a pixels center relative to the defined center
	   of this operation.
	 - output: the a,b,c, and d parameters for the points of intersection
	           it also checks that the ray points toward the pixel, 
	           if not, all variables are set to -1 
	 - attention: this routine can only be applied if the defined center 
	              for this operation is not surrounded by the four pixel edges!!!!
	
	 What are a,b,c, and d?
	 They represent the four different sides of a pixel and measure the distance 
	 from a corner to the point of intersection in units of pixel edges.
	 A value between 0 and 1 represents an intersection. The upper left corer 
	 is the startingpoint of a and then clockwise we get the startingpoints of
	 the values b, c, and d. The direction in which they measure distances is also 
	 clockwise:
	
			|-------a---->---
	       7               |
	       |     Pixel     |
	       |               b
	       d   w/ len=1    |
	       |   w/ area=1   |
	       |               v
			---<----c-------|
	"""
	
	# prepare dictionary; the I symbolizes the circle
	dI = {}
	
	# calculate the slope parameters
	one_div_m = np.tan(phi)

	if one_div_m!=0.:
		m = 1./one_div_m
		
		a = (ce[1]+0.5) * one_div_m - (ce[0]-0.5)
		b = - (ce[0]+0.5) * m + (ce[1]+0.5)
		c = - (ce[1]-0.5) * one_div_m + (ce[0]+0.5)
		d = (ce[0]-0.5) * m - (ce[1]-0.5)
		
	else:
		a = - (ce[0]-0.5)
		b = -1
		c = (ce[0]+0.5)
		d = -1
	
	# the I sign symbolises the ray
	# we check first if the number is on the edge (ok) and if so, we check, 
	# whether the solution is actually in the direction of the ray
	# and not in its negative direction
	if ok(a): 
		if ((ce[0]-0.5+a) * np.sin(phi) + (ce[1]+0.5) * np.cos(phi) >= 0.): dI["a_I"] = a
	if ok(b): 
		if ((ce[0]+0.5) * np.sin(phi) + (ce[1]+0.5-b) * np.cos(phi) >= 0.): dI["b_I"] = b
	if ok(c): 
		if ((ce[0]+0.5-c) * np.sin(phi) + (ce[1]-0.5) * np.cos(phi) >= 0.): dI["c_I"] = c
	if ok(d): 
		if ((ce[0]-0.5) * np.sin(phi) + (ce[1]-0.5+d) * np.cos(phi) >= 0.): dI["d_I"] = d
		
	return dI

def is_geq_than(x,y):
	""" x is not None and greater than or equal to y """
	if x != None: 
		if x>=y: return True
	
	return False

def is_leq_than(x,y):
	""" x is not None and less than or equal to y """
	if x != None: 
		if x<=y: return True
	
	return False
	
def is_neq_than(x,y):
	""" x is not None and not equal to y """
	if x != None: 
		if x!=y: return True
	
	return False
	

def get_area_summed(pos,c,phi,r0,do,dI,phi_c,R_corner,Phi_corner,X_rim,ring_pos):
	""" - given a pixel, a ray, and a circle including all its properties,
		  this function determines the intersection area of the pixel and
		  the area inside the radius of the circle and left to the ray
		- phi_c is the phi value correspondng to the center of the pixel
		- r_corner is an array containing the 4 radial distances to the 
		  corners of the pixel (r1,r2,r3,r4)
	"""
	
	if pos == "TR": return get_area_summed_TR(c,phi,r0,do,dI,phi_c,R_corner)
	elif pos == "T": return get_area_summed_T(c,phi,r0,do,dI,phi_c,R_corner)
	elif pos == "M": return get_area_summed_M(c,phi,r0,do,dI,phi_c,R_corner,Phi_corner,X_rim,ring_pos)
	else: 
		print("Error! Unexpected position for area determination:",pos)
		sys.exit()
	

def get_area_summed_TR(c,phi,r0,do,dI,phi_c,R_corner):
	""" given a pixel in the position top right ('TR'), a ray, 
		and a circle including all its properties,
		this function determines the intersection area of the pixel and
		the area inside the radius of the circle and left to the ray """
	
	if do.get("d_o") != None:
		if do.get("c_o") != None:
			# we have do and co

			if phi <= do.get("phi_d_o"): return 0.
			elif phi >= do.get("phi_c_o"): return 0.5 * ( (do.get("phi_c_o")-do.get("phi_d_o"))*r0*r0 - (c[0]-0.5)*do.get("d_o") - (c[1]-0.5)*(1.-do.get("c_o")) ) 
			elif is_leq_than(dI.get("d_I"),do.get("d_o")): return 0.5 * ( (phi-do.get("phi_d_o"))*r0*r0 - (c[0]-0.5)*(do.get("d_o")-dI.get("d_I")) )
			elif is_geq_than(dI.get("c_I"),do.get("c_o")): return 0.5 * ( (phi-do.get("phi_d_o"))*r0*r0 - (c[0]-0.5)*do.get("d_o") - (c[1]-0.5)*(1.-dI.get("c_I")) ) 
			else: return 0.5 * ( (do.get("phi_c_o")-do.get("phi_d_o"))*r0*r0 - (c[0]-0.5)*do.get("d_o") - (c[1]-0.5)*(1.-do.get("c_o")) ) 
						
			
		elif do.get("b_o") != None:
			# we have do and bo

			if phi <= do.get("phi_d_o"): return 0.
			elif is_leq_than(dI.get("d_I"),do.get("d_o")): return 0.5 * ( (phi-do.get("phi_d_o"))*r0*r0 - (c[0]-0.5)*(do.get("d_o")-dI.get("d_I")) )
			elif ((dI.get("c_I")!=None) and (phi<=do.get("phi_b_o"))): return 0.5 * ( (phi-do.get("phi_d_o"))*r0*r0 - (c[0]-0.5)*do.get("d_o") - (c[1]-0.5)*(1.-dI.get("c_I")) ) 
			elif dI.get("c_I")!=None: return 0.5 * ( (do.get("phi_b_o")-do.get("phi_d_o"))*r0*r0 - (c[0]-0.5)*do.get("d_o") - (c[1]-0.5)*(1.-dI.get("c_I")) + (c[0]+0.5)*(dI.get("b_I")-do.get("b_o")) ) 
			else: return 0.5 * ( (do.get("phi_b_o")-do.get("phi_d_o"))*r0*r0 - (c[0]-0.5)*do.get("d_o") - (c[1]-0.5) + (c[0]+0.5)*(1.-do.get("b_o")) ) 
			
		else: print("Error! 1. Unexpected case in TR area determination:")

			
	elif do.get("a_o") != None:
		if do.get("c_o") != None:
			# we have ao and co

			if ((not dI.get("d_I")) and (not dI.get("c_I")) and (phi <= do.get("phi_a_o"))): return 0.
			elif ((dI.get("d_I")!=None) and (phi<=do.get("phi_a_o"))): return 0.5 * (1.-dI.get("d_I")) * dI.get("a_I")
			elif dI.get("d_I")!=None: return 0.5 * ( (c[1]+0.5)*do.get("a_o") + (phi-do.get("phi_a_o"))*r0*r0 - (c[0]-0.5)*(1.-dI.get("d_I")) )
			elif is_geq_than(dI.get("c_I"),do.get("c_o")): return 0.5 * ( (c[1]+0.5)*do.get("a_o") + (phi-do.get("phi_a_o"))*r0*r0 - (c[0]-0.5) - (c[1]-0.5)*(1.-dI.get("c_I")) )
			else: return 0.5 * ( (c[1]+0.5)*do.get("a_o") + (do.get("phi_c_o")-do.get("phi_a_o"))*r0*r0 - (c[0]-0.5) - (c[1]-0.5)*(1.-do.get("c_o")) )
			
		elif do.get("b_o") != None:
			# we have ao and bo

			if ((not dI.get("d_I")) and (not dI.get("c_I")) and (phi <= do.get("phi_a_o"))): return 0.
			elif (dI.get("d_I")!=None):
				if (phi<=do.get("phi_a_o")): return 0.5 * (1.-dI.get("d_I")) * dI.get("a_I")
				elif (phi<=do.get("phi_b_o")): return 0.5 * ( (c[1]+0.5)*do.get("a_o") + (phi-do.get("phi_a_o"))*r0*r0 - (c[0]-0.5)*(1.-dI.get("d_I")) )
				elif (dI.get("b_I")!=None): return 0.5 * ( (c[1]+0.5)*do.get("a_o") + (do.get("phi_b_o")-do.get("phi_a_o"))*r0*r0 + (c[0]+0.5)*(dI.get("b_I")-do.get("b_o")) - (c[0]-0.5)*(1.-dI.get("d_I")) )
				elif (dI.get("b_I")!=None): return 0.5 * ( (c[1]+0.5)*do.get("a_o") + (do.get("phi_b_o")-do.get("phi_a_o"))*r0*r0 + (c[0]+0.5)*(dI.get("b_I")-do.get("b_o")) - (c[0]-0.5)*(1.-dI.get("d_I")) )
				elif dI.get("d_I")==0.: return 0.5 * ( (c[1]+0.5)*do.get("a_o") + (do.get("phi_b_o")-do.get("phi_a_o"))*r0*r0 - (c[0]-0.5) + (c[0]+0.5)*(1.-do.get("b_o")) )
				else: print("Error! 2.1 Unexpected case.")
			elif (dI.get("c_I")!=None):
				if (phi<=do.get("phi_a_o")): return 0.5 * (1. - dI.get("c_I") + dI.get("a_I"))
				elif (phi<=do.get("phi_b_o")): return 0.5 * ( (c[1]+0.5)*do.get("a_o") + (phi-do.get("phi_a_o"))*r0*r0 - (c[0]-0.5) - (c[1]-0.5)*(1.-dI.get("c_I")) )
				elif dI.get("b_I")!=None: return 0.5 * ( (c[1]+0.5)*do.get("a_o") + (do.get("phi_b_o")-do.get("phi_a_o"))*r0*r0 - (c[0]-0.5) - (c[1]-0.5)*(1.-dI.get("c_I")) + (c[0]+0.5)*(dI.get("b_I")-do.get("b_o")) )
				elif dI.get("c_I")==1.: return 0.5 * ( (c[1]+0.5)*do.get("a_o") + (do.get("phi_b_o")-do.get("phi_a_o"))*r0*r0 - (c[0]-0.5) + (c[0]+0.5)*(1.-do.get("b_o")) )
				else: print("Error! 2.2 Unexpected case.")
			else: return 0.5 * ( (c[1]+0.5)*do.get("a_o") + (do.get("phi_b_o")-do.get("phi_a_o"))*r0*r0 - (c[0]-0.5) - (c[1]-0.5) + (c[0]+0.5)*(1.-do.get("b_o")) )
		
		else: print("Error! 2. Unexpected case in TR area determination:")

	else:
		# circle does not intersect pixel
		if r0<=R_corner[2]: 
			# pixel is outside the circle

			return 0.
			
		elif r0>=R_corner[0]:
			# pixel is inside the circle

			if ((not dI.get("d_I")) and (not dI.get("c_I"))):
				if (phi <= phi_c): return 0.
				else: return 1.
			
			if (dI.get("d_I")!=None):		# here we want to catch the case that d=0 and c=1, which corresponds to the center being in the botton left corner
				if dI.get("a_I")!=None: return 0.5 * (1.-dI.get("d_I")) * dI.get("a_I")
				elif dI.get("b_I")!=None: return 0.5 * (1.-dI.get("d_I") + dI.get("b_I"))
				elif dI.get("d_I")==0.: # here we want to catch the case that d=0 and c=1, which corresponds to the center being in the botton left corner
					if phi <= phi_c: return 0.
					else: return 1.				
				else:print("Error! 3.1 Unexpected case in TR area determination:")
			elif dI.get("c_I")!=None:		# here we want to catch the case that d=0 and c=1, which corresponds to the center being in the botton left corner
				if dI.get("a_I")!=None: return 0.5 * (1.-dI.get("c_I") + dI.get("a_I"))
				elif dI.get("b_I")!=None: return 1. - 0.5 * (1.-dI.get("b_I")) * dI.get("c_I")
				elif dI.get("c_I")==1.: # here we want to catch the case that d=0 and c=1, which corresponds to the center being in the botton left corner
					if phi <= phi_c: return 0.
					else: return 1.	
				else:print("Error! 3.2 Unexpected case in TR area determination:")
			else: return 1.
			
		else: print("Error! 4. Unexpected case in TR area determination:")
		
	print("c,phi,r0,do,dI,phi_c,R_corner")
	print(c,phi,r0,do,dI,phi_c,R_corner)
	sys.exit()

def get_area_summed_T(c,phi,r0,do,dI,phi_c,R_corner):
	""" given a pixel in the position top ('T'), a ray, 
		and a circle including all its properties,
		this function determines the intersection area of the pixel and
		the area inside the radius of the circle and left to the ray"""

	
	if r0 <= c[1]-0.5: 
		# case pixel outside the radius
		return 0.
	
	elif do["has_c_o"]==2: 
		# case with two c solutions
		if pi<phi<=do["phi_c_o+"]: return 0.
		elif phi>do["phi_c_o+"]: return 0.5 * ( (phi-do["phi_c_o+"])*r0*r0 - (c[1]-0.5)*(do.get("c_o+")-dI.get("c_I")) )
		elif phi<=do["phi_c_o"]: return 0.5 * ( (pix2 + phi-do["phi_c_o+"])*r0*r0 - (c[1]-0.5)*(do.get("c_o+")-dI.get("c_I")) )
		else: return 0.5 * ( (pix2 + do["phi_c_o"] - do["phi_c_o+"])*r0*r0 - (c[1]-0.5)*(do.get("c_o+")-do.get("c_o")) )
		
	elif do["has_c_o"]==1:
		if do.get("d_o") != None:
			# case cd
			if (not dI.get("c_I")): 
				if (phi>pi): return 0.
				else: return 0.5 * ( (pix2 + do["phi_c_o"] - do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") - (1.-do.get("c_o"))*(c[1]-0.5) )
			elif is_leq_than(dI.get("d_I"),do.get("d_o")): return 0.5 * (1.-dI.get("c_I"))*dI.get("d_I")
			elif phi>=do.get("phi_d_o"): return 0.5 * ( (phi-do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") - (1.-dI.get("c_I"))*(c[1]-0.5) )
			elif phi<=do.get("phi_c_o"): return 0.5 * ( (pix2 + phi-do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") - (1.-dI.get("c_I"))*(c[1]-0.5) )
			else: return 0.5 * ( (pix2 + do["phi_c_o"] - do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") - (1.-do.get("c_o"))*(c[1]-0.5) )
			
		elif do.get("b_o") != None:
			# case cb
			if pi<phi<=do["phi_c_o"]: return 0.
			elif phi>=do.get("phi_c_o"): return 0.5 * ( (phi-do.get("phi_c_o"))*r0*r0 - (c[1]-0.5)*(do.get("c_o")-dI.get("c_I")) )
			elif phi<=do.get("phi_b_o"): return 0.5 * ( (pix2 + phi-do.get("phi_c_o"))*r0*r0 - (c[1]-0.5)*(do.get("c_o")-dI.get("c_I")) )
			elif is_geq_than(dI.get("b_I"),do.get("b_o")): return 0.5 * ( (pix2 + do.get("phi_b_o")-do.get("phi_c_o"))*r0*r0 + (c[0]+0.5)*(dI.get("b_I")-do.get("b_o")) - (c[1]-0.5)*(do.get("c_o")-dI.get("c_I")) )
			else: return 0.5 * ( (pix2 + do.get("phi_b_o")-do.get("phi_c_o"))*r0*r0 + (c[0]+0.5)*(1.-do.get("b_o")) - (c[1]-0.5)*do.get("c_o") )
			
		else: print("Error! Here: get_area_summed_T; has_c_o == 1")
			
	elif ((do.get("b_o") != None) and (do.get("d_o") != None) and (not do.get("a_o"))):
		# case bd and not a (practivally)
		if (not dI.get("c_I")): 
			if (phi>pi): return 0.
			else: return 0.5 * ( (pix2 + do["phi_b_o"] - do["phi_d_o"])*r0*r0 + (0.5+c[0])*(1.-do.get("b_o")) + (0.5-c[0])*do.get("d_o") - (c[1]-0.5) )
		elif is_leq_than(dI.get("d_I"),do.get("d_o")): return 0.5 * (1.-dI.get("c_I"))*dI.get("d_I")
		elif phi>=do.get("phi_d_o"): return 0.5 * ( (phi-do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") - (c[1]-0.5)*(1.-dI.get("c_I")) )
		elif phi<=do.get("phi_b_o"): return 0.5 * ( (pix2 + phi-do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") - (c[1]-0.5)*(1.-dI.get("c_I")) )
		elif is_geq_than(dI.get("b_I"),do.get("b_o")): return 0.5 * ( (pix2 + do["phi_b_o"]-do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") + (0.5+c[0])*(dI.get("b_I")-do.get("b_o")) - (c[1]-0.5)*(1.-dI.get("c_I")) )
		else: return 0.5 * ( (pix2 + do["phi_b_o"] - do["phi_d_o"])*r0*r0 + (0.5+c[0])*(1.-do.get("b_o")) + (0.5-c[0])*do.get("d_o") - (c[1]-0.5) )
					
	elif do["has_a_o"]==2:
		# case two a solutions (this should imply that b and d also have solutions)

		if (not dI.get("c_I")): 
			if (phi>pi): return 0.
			else: return 0.5 * ( (do["phi_a_o"]-do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") + (c[1]+0.5)*(do.get("a_o+")-do.get("a_o")) + (do["phi_b_o"]-do["phi_a_o+"])*r0*r0 + (0.5+c[0])*(1.-do.get("b_o")) - (c[1]-0.5) )
		elif is_leq_than(dI.get("d_I"),do.get("d_o")): return 0.5 * (1.-dI.get("c_I"))*dI.get("d_I")
		elif do.get("phi_d_o")<=phi<=do.get("phi_a_o"): return 0.5 * ( (phi-do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") - (c[1]-0.5)*(1.-dI.get("c_I")) )
		elif ((do.get("phi_a_o")<=phi) or (phi<=do.get("phi_a_o+"))): 
			# ~ print(">>>",phi,do["phi_a_o"],do["phi_a_o+"],"asd",(do.get("phi_a_o")<=phi),(phi<=do.get("phi_a_o+")),"asdasd",do["phi_d_o"],do.get("d_o"),dI.get("a_I"),do.get("a_o"),dI.get("c_I"))
			return 0.5 * ( (do["phi_a_o"]-do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") + (c[1]+0.5)*(dI.get("a_I")-do.get("a_o")) - (c[1]-0.5)*(1.-dI.get("c_I")) )
		elif do.get("phi_a_o+")<=phi<=do.get("phi_b_o"): return 0.5 * ( (do["phi_a_o"]-do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") + (c[1]+0.5)*(do.get("a_o+")-do.get("a_o")) + (phi-do["phi_a_o+"])*r0*r0 - (c[1]-0.5)*(1.-dI.get("c_I")) )
		elif is_geq_than(dI.get("b_I"),do.get("b_o")): return 0.5 * ( (do["phi_a_o"]-do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") + (c[1]+0.5)*(do.get("a_o+")-do.get("a_o")) + (do["phi_b_o"]-do["phi_a_o+"])*r0*r0 + (0.5+c[0])*(dI.get("b_I")-do.get("b_o")) - (c[1]-0.5)*(1.-dI.get("c_I")) )
		else: return 0.5 * ( (do["phi_a_o"]-do["phi_d_o"])*r0*r0 + (0.5-c[0])*do.get("d_o") + (c[1]+0.5)*(do.get("a_o+")-do.get("a_o")) + (do["phi_b_o"]-do["phi_a_o+"])*r0*r0 + (0.5+c[0])*(1.-do.get("b_o")) - (c[1]-0.5) )
		
	elif do["has_a_o"]==1:
		if do.get("b_o") != None:
			# case ab, this case should imply that there is no d sollution

			if (not dI.get("c_I")): 
				if (phi>pi): return 0.
				else: return 0.5 * ( -(c[1]-0.5) + (0.5-c[0]) + (c[1]+0.5)*do.get("a_o") + (do["phi_b_o"]-do["phi_a_o"])*r0*r0 + (0.5+c[0])*(1.-do.get("b_o")) )
			elif dI.get("d_I")!=None: return 0.5 * (1.-dI.get("c_I"))*dI.get("d_I")
			elif is_leq_than(dI.get("a_I"),do.get("a_o")): return 0.5 * ( -(c[1]-0.5)*(1.-dI.get("c_I")) + (0.5-c[0]) + (c[1]+0.5)*dI.get("a_I") )
			elif do.get("phi_a_o")<=phi<=do.get("phi_b_o"): return 0.5 * ( -(c[1]-0.5)*(1.-dI.get("c_I")) + (0.5-c[0]) + (c[1]+0.5)*do.get("a_o") + (phi-do["phi_a_o"])*r0*r0 )
			elif is_geq_than(dI.get("b_I"),do.get("b_o")): return 0.5 * ( -(c[1]-0.5)*(1.-dI.get("c_I")) + (0.5-c[0]) + (c[1]+0.5)*do.get("a_o") + (do["phi_b_o"]-do["phi_a_o"])*r0*r0 + (0.5+c[0])*(dI.get("b_I")-do.get("b_o")) )
			else: return 0.5 * ( -(c[1]-0.5) + (0.5-c[0]) + (c[1]+0.5)*do.get("a_o") + (do["phi_b_o"]-do["phi_a_o"])*r0*r0 + (0.5+c[0])*(1.-do.get("b_o")) )
			
			
		
		elif do.get("d_o") != None:
			# case ad, this case should imply that there is no b sollution

			if (not dI.get("c_I")): 
				if (phi>pi): return 0.
				else: return 0.5 * ( -(c[1]-0.5) + (0.5-c[0])*do.get("d_o") + (do.get("phi_a_o")-do.get("phi_d_o"))*r0*r0 + (c[1]+0.5)*(1.-do.get("a_o")) + (0.5+c[0]) )
			elif is_leq_than(dI.get("d_I"),do.get("d_o")): return 0.5 * (1.-dI.get("c_I"))*dI.get("d_I")
			elif do.get("phi_d_o")<=phi<=do.get("phi_a_o"): return 0.5 * ( -(c[1]-0.5)*(1.-dI.get("c_I")) + (0.5-c[0])*do.get("d_o") + (phi-do.get("phi_d_o"))*r0*r0 )
			elif is_geq_than(dI.get("a_I"),do.get("a_o")): return 0.5 * ( -(c[1]-0.5)*(1.-dI.get("c_I")) + (0.5-c[0])*do.get("d_o") + (do.get("phi_a_o")-do.get("phi_d_o"))*r0*r0 + (c[1]+0.5)*(dI.get("a_I")-do.get("a_o")) )
			elif dI.get("b_I")!=None: return 0.5 * ( -(c[1]-0.5)*(1.-dI.get("c_I")) + (0.5-c[0])*do.get("d_o") + (do.get("phi_a_o")-do.get("phi_d_o"))*r0*r0 + (c[1]+0.5)*(1.-do.get("a_o")) + (0.5+c[0])*dI.get("b_I") )
			else: return 0.5 * ( -(c[1]-0.5) + (0.5-c[0])*do.get("d_o") + (do.get("phi_a_o")-do.get("phi_d_o"))*r0*r0 + (c[1]+0.5)*(1.-do.get("a_o")) + (0.5+c[0]) )
			
			
		else: print("Error! Here: get_area_summed_T; has_a_o == 1")		
	
	elif ((r0 >= R_corner[0]) and (r0 >= R_corner[3])):
		# case pixel fully inside the radius

		if (not dI.get("c_I")): 
			if (phi>pi): return 0.
			else: return 1. 
		elif dI.get("d_I")!=None: return 0.5 * (1.-dI.get("c_I"))*dI.get("d_I")
		elif dI.get("a_I")!=None: return 0.5 * (1.-dI.get("c_I")+dI.get("a_I"))
		elif dI.get("b_I")!=None: return 1. - 0.5 * (1.-dI.get("b_I"))*dI.get("c_I")
		else: 
			# if we are sitting on an edge this case is relevant
			if (phi>pi): return 0.
			else: return 1.

	
	else: print("Error! End of get_area_summed_T, couldn't find a suitable case.")
	
	print("c,phi,r0,do,dI,phi_c,R_corner")
	print(c,phi,r0,do,dI,phi_c,R_corner)
	sys.exit()


def get_area_summed_M(c,phi,r0,do,dI,phi_c,R_corner,Phi_corner,X_rim,ring_pos):
	""" given a pixel in the position middle ('M'), a ray, 
		and a circle including all its properties,
		this function determines the intersection area of the pixel and
		the area inside the radius of the circle and left to the ray;
		this function, oposed to the others, is not derived on paper. 
		I was to impatient and just wrote it down in python. """
	
	# the area will be increased during this function
	area = 0
	
	# we are going in phi direction over the pixel area and add up its contributions
	# until we reach the phi value
	idx = 0
	if phi<=Phi_corner[idx]:
		if ring_pos[idx] == "p": area += 0.5 * ( (dI["a_I"] - X_rim[(idx-1)%4])*(0.5+c[1]) )
		elif ring_pos[idx] == "c": area += 0.5 * ( phi*r0*r0 )
		elif ring_pos[idx] == "cp": area += 0.5 * ( (dI["a_I"] - X_rim[(idx-1)%4])*(0.5+c[1]) )
		elif ring_pos[idx] == "pc": 
			if phi<do["phi_a_o"]: area += 0.5 * ( (dI["a_I"] - X_rim[(idx-1)%4])*(0.5+c[1]) )
			else: area += 0.5 * ( (do["a_o"] - X_rim[(idx-1)%4])*(0.5+c[1]) + (phi-do["phi_a_o"])*r0*r0 )
		elif ring_pos[idx] == "cpc": 
			if phi<do["phi_a_o+"]: area += 0.5 * ( (dI["a_I"] - X_rim[(idx-1)%4])*(0.5+c[1]) )
			else: area += 0.5 * ( (do["a_o+"] - X_rim[(idx-1)%4])*(0.5+c[1]) + (phi-do["phi_a_o+"])*r0*r0 )
		return area
	else:
		if ring_pos[idx] == "p": area += 0.5 * ( (1. - X_rim[(idx-1)%4])*(0.5+c[1]) )
		elif ring_pos[idx] == "c": area += 0.5 * ( Phi_corner[idx]*r0*r0)
		elif ring_pos[idx] == "cp": area += 0.5 * ( (1. - X_rim[(idx-1)%4])*(0.5+c[1]) )
		elif ring_pos[idx] == "pc": area += 0.5 * ( (do["a_o"] - X_rim[(idx-1)%4])*(0.5+c[1]) + (Phi_corner[idx]-do["phi_a_o"])*r0*r0 )
		elif ring_pos[idx] == "cpc": area += 0.5 * ( (do["a_o+"] - X_rim[(idx-1)%4])*(0.5+c[1]) + (Phi_corner[idx]-do["phi_a_o+"])*r0*r0 )

	idx = 1
	if phi<=Phi_corner[idx]:
		if ring_pos[idx] == "p": area += 0.5 * ( dI["b_I"]*(0.5+c[0]) )
		elif ring_pos[idx] == "c": area += 0.5 * ( (phi-Phi_corner[idx-1])*r0*r0 )
		elif ring_pos[idx] == "cp": 
			if phi<do["phi_b_o"]: area += 0.5 * ( (phi-Phi_corner[idx-1])*r0*r0 )
			else: area += 0.5 * ( (do["phi_b_o"]-Phi_corner[idx-1])*r0*r0 + (dI["b_I"] - do["b_o"])*(0.5+c[0]) )
		elif ring_pos[idx] == "pc": 
			if phi<do["phi_b_o"]: area += 0.5 * ( dI["b_I"]*(0.5+c[0]) )
			else: area += 0.5 *  ( do["b_o"]*(0.5+c[0]) + (phi-do["phi_b_o"])*r0*r0 )
		elif ring_pos[idx] == "cpc": 
			if phi<do["phi_b_o"]: area += 0.5 * ( (phi-Phi_corner[idx-1])*r0*r0 )
			elif phi<do["phi_b_o+"]: area += 0.5 * ( (do["phi_b_o"]-Phi_corner[idx-1])*r0*r0 + (dI["b_I"] - do["b_o"])*(0.5+c[0]) )
			else: area += 0.5 * ( (do["phi_b_o"]-Phi_corner[idx-1])*r0*r0 + (do["b_o+"] - do["b_o"])*(0.5+c[0]) + (phi-do["phi_b_o+"])*r0*r0 )
		return area
	else:
		if ring_pos[idx] == "p": area += 0.5 * ( (0.5+c[0]) )
		elif ring_pos[idx] == "c": area += 0.5 * ( (Phi_corner[idx]-Phi_corner[idx-1]) * r0*r0)
		elif ring_pos[idx] == "cp": area += 0.5 * ( (do["phi_b_o"]-Phi_corner[idx-1])*r0*r0 + (1. - do["b_o"])*(0.5+c[0]) )
		elif ring_pos[idx] == "pc": area += 0.5 * ( do["b_o"]*(0.5+c[0]) + (Phi_corner[idx]-do["phi_b_o"])*r0*r0 )
		elif ring_pos[idx] == "cpc": area += 0.5 * ( (do["phi_b_o"]-Phi_corner[idx-1])*r0*r0 + (do["b_o+"] - do["b_o"])*(0.5+c[0]) + (Phi_corner[idx]-do["phi_b_o+"])*r0*r0 )

	idx = 2
	if phi<=Phi_corner[idx]:
		if ring_pos[idx] == "p": area += 0.5 * ( dI["c_I"]*(0.5-c[1]) )
		elif ring_pos[idx] == "c": area += 0.5 * ( (phi-Phi_corner[idx-1])*r0*r0 )
		elif ring_pos[idx] == "cp": 
			if phi<do["phi_c_o"]: area += 0.5 * ( (phi-Phi_corner[idx-1])*r0*r0 )
			else: area += 0.5 * ( (do["phi_c_o"]-Phi_corner[idx-1])*r0*r0 + (dI["c_I"] - do["c_o"])*(0.5-c[1]) )
		elif ring_pos[idx] == "pc": 
			if phi<do["phi_c_o"]: area += 0.5 * ( dI["c_I"]*(0.5-c[1]) )
			else: area += 0.5 *  ( do["c_o"]*(0.5-c[1]) + (phi-do["phi_c_o"])*r0*r0 )
		elif ring_pos[idx] == "cpc": 
			if phi<do["phi_c_o"]: area += 0.5 * ( (phi-Phi_corner[idx-1])*r0*r0 )
			elif phi<do["phi_c_o+"]: area += 0.5 * ( (do["phi_c_o"]-Phi_corner[idx-1])*r0*r0 + (dI["c_I"] - do["c_o"])*(0.5-c[1]) )
			else: area += 0.5 * ( (do["phi_c_o"]-Phi_corner[idx-1])*r0*r0 + (do["c_o+"] - do["c_o"])*(0.5-c[1]) + (phi-do["phi_c_o+"])*r0*r0 )
		return area
	else:
		if ring_pos[idx] == "p": area += 0.5 * ( (0.5-c[1]) )
		elif ring_pos[idx] == "c": area += 0.5 * ( (Phi_corner[idx]-Phi_corner[idx-1]) * r0*r0)
		elif ring_pos[idx] == "cp": area += 0.5 * ( (do["phi_c_o"]-Phi_corner[idx-1])*r0*r0 + (1. - do["c_o"])*(0.5-c[1]) )
		elif ring_pos[idx] == "pc": area += 0.5 * ( do["c_o"]*(0.5-c[1]) + (Phi_corner[idx]-do["phi_c_o"])*r0*r0 )
		elif ring_pos[idx] == "cpc": area += 0.5 * ( (do["phi_c_o"]-Phi_corner[idx-1])*r0*r0 + (do["c_o+"] - do["c_o"])*(0.5-c[1]) + (Phi_corner[idx]-do["phi_c_o+"])*r0*r0 )

	idx = 3
	if phi<=Phi_corner[idx]:
		if ring_pos[idx] == "p": area += 0.5 * ( dI["d_I"]*(0.5-c[0]) )
		elif ring_pos[idx] == "c": area += 0.5 * ( (phi-Phi_corner[idx-1])*r0*r0 )
		elif ring_pos[idx] == "cp": 
			if phi<do["phi_d_o"]: area += 0.5 * ( (phi-Phi_corner[idx-1])*r0*r0 )
			else: area += 0.5 * ( (do["phi_d_o"]-Phi_corner[idx-1])*r0*r0 + (dI["d_I"] - do["d_o"])*(0.5-c[0]) )
		elif ring_pos[idx] == "pc": 
			if phi<do["phi_d_o"]: area += 0.5 * ( dI["d_I"]*(0.5-c[0]) )
			else: area += 0.5 *  ( do["d_o"]*(0.5-c[0]) + (phi-do["phi_d_o"])*r0*r0 )
		elif ring_pos[idx] == "cpc": 
			if phi<do["phi_d_o"]: area += 0.5 * ( (phi-Phi_corner[idx-1])*r0*r0 )
			elif phi<do["phi_d_o+"]: area += 0.5 * ( (do["phi_d_o"]-Phi_corner[idx-1])*r0*r0 + (dI["d_I"] - do["d_o"])*(0.5-c[0]) )
			else: area += 0.5 * ( (do["phi_d_o"]-Phi_corner[idx-1])*r0*r0 + (do["d_o+"] - do["d_o"])*(0.5-c[0]) + (phi-do["phi_d_o+"])*r0*r0 )
		return area
	else:
		if ring_pos[idx] == "p": area += 0.5 * ( (0.5-c[0]) )
		elif ring_pos[idx] == "c": area += 0.5 * ( (Phi_corner[idx]-Phi_corner[idx-1]) * r0*r0)
		elif ring_pos[idx] == "cp": area += 0.5 * ( (do["phi_d_o"]-Phi_corner[idx-1])*r0*r0 + (1. - do["d_o"])*(0.5-c[0]) )
		elif ring_pos[idx] == "pc": area += 0.5 * ( do["d_o"]*(0.5-c[0]) + (Phi_corner[idx]-do["phi_d_o"])*r0*r0 )
		elif ring_pos[idx] == "cpc": area += 0.5 * ( (do["phi_d_o"]-Phi_corner[idx-1])*r0*r0 + (do["d_o+"] - do["d_o"])*(0.5-c[0]) + (Phi_corner[idx]-do["phi_d_o+"])*r0*r0 )
	
	
	idx = 0
	if ring_pos[idx] == "p": area += 0.5 * ( dI["a_I"]*(0.5+c[1]) )
	elif ring_pos[idx] == "c": area += 0.5 * ( (phi-Phi_corner[(idx-1)%4])*r0*r0 )
	elif ring_pos[idx] == "cp": 
		if phi<do["phi_a_o"]: area += 0.5 * ( (phi-Phi_corner[(idx-1)%4])*r0*r0 )
		else: area += 0.5 * ( (do["phi_a_o"]-Phi_corner[(idx-1)%4])*r0*r0 + (dI["a_I"] - do["a_o"])*(0.5+c[1]) )
	elif ring_pos[idx] == "pc": area += 0.5 * ( dI["a_I"]*(0.5+c[1]) )
	elif ring_pos[idx] == "cpc": 
		if phi<do["phi_a_o"]: area += 0.5 * ( (phi-Phi_corner[(idx-1)%4])*r0*r0 )
		else: area += 0.5 * ( (do["phi_a_o"]-Phi_corner[(idx-1)%4])*r0*r0 + (dI["a_I"] - do["a_o"])*(0.5+c[1]) )
	return area



def summed_to_indivudual_contribution(sumd):
	""" given a summed contribution array, it gives the indivudual contributions """
	return (sumd - np.pad(sumd,((0,0),(1,0)), mode='constant')[:, :-1] - np.pad(sumd,((1,0),(0,0)), mode='constant')[:-1, :] + np.pad(sumd,((1,0),(1,0)), mode='constant')[:-1, :-1])

def get_rotated_pixel_contributionn(c,pos,n_phi_min,n_phi_max,n_r_min, n_r_max, N_phi,N_r, dphi, dr):
	""" given a pixel's position and the new grid's pixels properties,
		this function calculates the contribution to all considered new grid cells.
		in order to work properly, the new grid cells have to be chosen appropriately.
	"""
	####################################################################
	# 1. generate output array that covers the complete cartesian pixel + TEST FOR DEBUGGING?
	####################################################################

	dn_phi = int(n_phi_max - n_phi_min + 1) if (n_phi_max >= n_phi_min) else int(n_phi_max + N_phi - n_phi_min + 1)
	dn_r = int(n_r_max - n_r_min + 1)
	contr = np.zeros((dn_phi,dn_r))
	debug = False

	####################################################################
	# 2. precalculate some variable
	####################################################################
	
	phi_c = np.arctan2(c[0],c[1])				# phi-value of the cart. pixel's center
	if phi_c < 0.: phi_c += pix2
	
	R_corner = np.zeros(4)						# vector containing distances to the pixel's corners (r1,r2,r3,r4)
	R_corner[0] = get_r(c + np.array([0.5,0.5]))
	R_corner[1] = get_r(c + np.array([0.5,-0.5]))
	R_corner[2] = get_r(c + np.array([-0.5,-0.5]))
	R_corner[3] = get_r(c + np.array([-0.5,0.5]))
	
	if n_phi_min > n_phi_max: 
		if pos != "T":
			print("Error! Unexpected position, given the contributing bin information: ",pos,n_phi_min,n_phi_max)
			sys.exit()
		phi_vals = [dphi*(n_phi_min+1+_) for _ in range(int(N_phi - 1 - n_phi_min))] + [dphi*_ for _ in range(int(n_phi_max + 2))]
	elif n_phi_max > N_phi -1: 
		print("Error! Unexpected bin number encountered: ",n_phi_max)
		sys.exit()
	elif n_phi_max == N_phi - 1: 
		if n_phi_max>n_phi_min: phi_vals = [dphi*(n_phi_min+1+_) for _ in range(dn_phi-1)] + [0.]
		else: phi_vals = [0.]
	else: 
		# the "good" case, i.e., n_phi_min <= n_phi_max < N_phi-1
		phi_vals = [dphi*(n_phi_min+1+_) for _ in range(dn_phi)]
	
	r0_vals = [dr*(n_r_min+1+_) for _ in range(dn_r)]
	
	# dictionaries with info about intersection of circle with pixel
	dos = []
	Ring_pos = []
	for r0 in r0_vals: dos.append(get_circle_pixel_intersect(c,r0))
	
	if pos == "M": 
		# calculate some additional information for the central pixel
		# x_... are distances to the pixel rim and [phi_0,phi_1,phi_2,phi_3]
		# are the angles to the corners. Ring_pos [[..,..,..,..] ... ] stores for all cicles 
		# and for all pixel edges the current circle position. 		
		# in the end it tells us (in phi direction) which thing is 
		# leading to the limit of the area, circle or pixel
		# circle inside pixel: "c"
		# circle outside pixel: "p"
		# first circle, then pixel: "cp"
		# first pixel, then circle: "pc"
		# circle, pixel, circle: "cpc"
		phi_vals = [dphi*(1+_) for _ in range(dn_phi)]
		
		Phi_corner = np.zeros(4)
		Phi_corner[0] = get_phi(c + np.array([0.5,0.5]))
		Phi_corner[1] = get_phi(c + np.array([0.5,-0.5]))
		Phi_corner[2] = get_phi(c + np.array([-0.5,-0.5]))
		Phi_corner[3] = get_phi(c + np.array([-0.5,0.5]))
		
		X_rim = np.array([0.5+c[1],0.5+c[0],0.5-c[1],0.5-c[0]])
		if not np.all(X_rim>0.): 
			print("Error! Unexpected position of center in Pixel:",X_rim)
			sys.exit
		for i,r0 in enumerate(r0_vals): Ring_pos.append(get_ring_pos(r0,c,R_corner,dos[i]))
	else:
		Phi_corner = None
		X_rim = None
		Ring_pos = [None]*dn_r
	
	####################################################################
	# 3. loop over pixels of the rotated slice of the new grid using indices i,j
	####################################################################
	

	for i,phi in enumerate(phi_vals):
		# get intersection information fpr the phi value and the pixel
		if pos == "TR": dI = get_ray_pixel_intersect_non_M(c,phi)
		elif pos == "T": dI = get_ray_pixel_intersect_non_M(c,phi)
		elif pos == "M":  dI = get_ray_pixel_intersect_M(c,phi,Phi_corner)
		else: 
			print("Error! Unexpected position in get_rotated_pixel_contribution: ",pos)
			sys.exit()
			
		for j,(r0,do,ring_pos) in enumerate(zip(r0_vals,dos,Ring_pos)):
			# now get the area
			sub_area = get_area_summed(pos,c,phi,r0,do,dI,phi_c,R_corner,Phi_corner,X_rim,ring_pos)
			contr[i,j] += sub_area

	if debug: contr_old = np.copy(contr)
	contr = summed_to_indivudual_contribution(contr)

	if debug:
		if ((not np.all(contr>-1.e-12)) or (not np.all(contr<1.+1.e-12))):
			print("error: ",pos)
			print(contr_old)
			print(contr)
			sys.exit()
			# ~ contr = np.zeros((dn_phi,dn_r))
		
		if abs(contr.sum()-1.)>1e-12:

			print("err")
			print(contr.sum())
			print(contr)
			sys.exit()

	return dn_phi,dn_r,contr
			
	
	
def read_input(path_to_input_file):
	""" Reads the input.dat file that is part of the output of the code when converting and storing the detector. """
	
	inp = {}
	
	with open(path_to_input_file,"r") as f:
		for i,line in enumerate(f):
			if i>1: 
				line = line.replace("\n", "")
				[key,val] = line.split("=")
				key = key.split()[0]

				if key[0]=="N": inp[key] = int(val)
				elif key[:6]=="center": 
					val = val.replace("[","")
					val = val.replace("]","")
					val = val.replace(",","")
					inp[key] = [float(s) for s in val.split()]
				elif key in ["dr","dphi","R_px_max"]: inp[key] = float(val)
				elif key == "progress_bar": inp[key] = bool(val)
				else: inp[key] = val
				
	return inp


def convert(data_cart_path, data_cart,N_phi,N_r,center_shift = [0.,0.], results_folder_path = None, save_name = None, R_px_max_in = None, progress_bar = False, overwrite_results = False):
	""" 
		Task: Convert a cartesian grid into a polar grid
		Method: In order to find the polar representation of an originally cartesian detector grid, this code overlays the detector with a chosen polar grid and calculates the polar pixel values. 
				To that end, every cartesian pixel is assumed to represent a homogeneous flux distribution across its full area. The associated value of a polar pixel is then the sum of all values 
				of overlapping cartsian pixels each weighted with its geometrical intersection area with the polar pixel. 
		Parameters:
		- input grid data_cart_path[x,y] given in unit Jy
		- output grid output[phi,r] is given in Jy
		- conserves overall flux
		- center_shift is givel in [x_px,y_px], where x_px and y_px are real numbers
		  that indicate how much the center is shifted from the physical center of the data array
		- N_phi and N_r are the (integer) numbers of angular and radial "cells" used for the polar detector, respectively
		- N_phi has to satisfy: N_phi%4 == 0 and N_phi > 0
		- N_r has to satisfy: N_r > 0
		- the phi angle is sampled linearly from 0 to 2pi
		- the radius by default samples linearly from 0 to 0.5*sqrt(N_x**2+N_y**2), where (N_x,N_y) = data_cart.shape
		- results_folder_path is the path describing where the output (data_polar) shall be saved in a separate folder 
		   with the name save_name. If None, it is not saved. 
		- R_px_max is the maximum radius up to which this is done; if None, it uses the whole data region
		- for data_cart_path[x,y] the angle phi=0 is defined in pos. y-direction and phi=pi/2. in pos. x-direction
		- overwrite_results allows to overwrite results for the same project. If False it will instead change the project name a bit and save the new results withouth deleting the former results.
	"""
	# check conditions for number of cells
	if ((N_phi%4 != 0) or (N_r <= 0) or (N_phi <= 0)): 
		print("N_phi has to satisfy: N_phi%4 == 0\nN_r has to satisfy: N_r > 0\nN_phi, N_r = ",N_phi, N_r)
		sys.exit()
		
	if results_folder_path!=None: 
		if not os.path.exists(results_folder_path): 
			print("Error! Results can only be saved in an already existing folder: ",results_folder_path)
			sys.exit()
	
	# get image size
	(N_x,N_y) = data_cart.shape
	
	# determine the center and R_px_max if needed
	center = np.array([N_x/2.,N_y/2.]) + np.array(center_shift)
	if R_px_max_in == None: R_px_max = get_R_px_max(center,N_x,N_y)
	else: R_px_max = R_px_max_in

	# calculate delta_phi and delta_r
	dphi = pix2/N_phi
	dr   = R_px_max/N_r
	
	# prepare output array
	data_polar = np.zeros((N_phi,N_r))
	
	if progress_bar: ii, ii_print, ii_tot = 0, max(int(N_x*N_y*0.01),1), N_x*N_y
	for i in range(N_x):
		for j in range(N_y):
			# update progress
			if progress_bar:
				ii+=1
				if ii%ii_print == 0: print(int(ii/ii_tot*100),"% done")
				sys.stdout.flush()

			# 1. get original pixel center position in new polar grid in cartesian coordinates
			c_px = np.array([i+0.5,j+0.5]) - center
			
			# 2. get overall position of pixel in grid, i.e., Top, Bottom, Left, Right, Middle
			pos = get_overall_pixel_position(c_px)
			
			# 3. find the contributing phi and r bins. 
			# 	 bin furthest to the "left/center" gets "min"-label; bin furthest to the "right/infinity" gets "max"-label;
			#	 note that the "min"-labeled value can exceed the "max"-labeled value if the region contains the pos. y-axis
			(n_phi_min, n_phi_max, n_r_min, n_r_max) = get_contr_bins(c_px,pos,N_phi,N_r,dphi,dr)
			# ~ print("old",c_px,n_phi_min, n_phi_max, n_r_min, n_r_max,pos)
			if  n_r_min==-1: continue
			
			# 4. rotate c_px as well as the contributing bins, s.t. the only possible
			#    rotated states are in ["M","T","TR"]; "M" is also rotated!
			(c_new,pos_new,n_phi_min_new,n_phi_max_new) = rotate_pixel(c_px,pos,n_phi_min,n_phi_max,N_phi)
			# ~ print("new",c_new,pos_new,n_phi_min_new,n_phi_max_new)
			
			# 5. get contribution from rotated pixel 
			dn_phi, dn_r, contr_rot = get_rotated_pixel_contributionn(c_new,pos_new,n_phi_min_new,n_phi_max_new,n_r_min, n_r_max, N_phi,N_r, dphi, dr)
			
			# 6. rotate back and add contribution 
			if n_phi_max>=n_phi_min: 
				data_polar[n_phi_min:n_phi_max+1,n_r_min:n_r_max+1] += data_cart[i,j]*contr_rot
			else:
				ddn_phi = int(N_phi-n_phi_min)
				data_polar[n_phi_min:,n_r_min:n_r_max+1] += data_cart[i,j]*contr_rot[:ddn_phi,:]
				data_polar[:dn_phi-ddn_phi,n_r_min:n_r_max+1] += data_cart[i,j]*contr_rot[ddn_phi:,:]

	

	if results_folder_path != None:
		# save simulation input and output array 
		
		ii = 1
		if save_name is None: save_name_fin = "test"
		else: save_name_fin = save_name
		
		if not overwrite_results:
			while os.path.exists(results_folder_path+"/"+save_name_fin):
				print("Warning! The requested save_name is already in use:",save_name)
				save_name_fin = save_name + "(" + str(ii) + ")"
				print("To not lose anything, save_name is set to:",save_name_fin)
				ii += 1
			os.mkdir(results_folder_path+"/"+save_name_fin)
		else:
			if not os.path.exists(results_folder_path+"/"+save_name_fin): os.mkdir(results_folder_path+"/"+save_name_fin)
		
		# saving intput
		with open(results_folder_path+"/"+save_name_fin+"/input.dat",'w') as f:
			f.write("Output data is accessed via np.load('path/to/output.npz')['polar']\n")
			f.write("------------------------------------------------------------------\n")
			f.write(f"data_cart_path = {str(data_cart_path)}\n")
			f.write(f"N_phi          = {str(N_phi)}\n")
			f.write(f"N_r            = {str(N_r)}\n")
			f.write(f"center_shift   = {str(center_shift)}\n")
			f.write(f"save_name      = {str(save_name)}\n")
			f.write(f"save_name_fin  = {str(save_name_fin)}\n")
			f.write(f"R_px_max       = {str(R_px_max)}\n")
			f.write(f"R_px_max_in    = {str(R_px_max_in)}\n")
			f.write(f"N_x            = {str(N_x)}\n")
			f.write(f"N_y            = {str(N_y)}\n")
			f.write(f"center         = {str(center)}\n")
			f.write(f"dphi           = {str(dphi)}\n")
			f.write(f"dr             = {str(dr)}\n")
			f.write(f"progress_bar   = {str(progress_bar)}\n")
		
		# saving coordinates
		with open(results_folder_path+"/"+save_name_fin+"/phi_border_values.dat",'w') as f:
			for phi in range(N_phi+1): f.write(f"{phi*dphi}\n")
		
		with open(results_folder_path+"/"+save_name_fin+"/r_border_values.dat",'w') as f:
			for r in range(N_r+1): f.write(f"{r*dr}\n")
		

		# saving output; Note: This can be loaded with np.load("path/to/output.npz")['polar']
		np.savez_compressed(results_folder_path+"/"+save_name_fin+"/output.npz",polar = data_polar)
		# ~ np.save(results_folder_path+"/"+save_name_fin+"/output.npy",data_polar)
		
		
		
	
	return data_polar




if __name__ == "__main__":
	
	print("This is a package to convert a cartesian detector grid into a polar detector grid taking into account the various possible cross section areas of cartesian and polar pixels.\n")
	print("In order to use this package, please check out the short description and the few examples provided on my github page. Don't worry, it is not difficult to use.")
	
	
	
	
	
	
	
	
	
