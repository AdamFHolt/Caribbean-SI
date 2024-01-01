#!/usr/bin/env python 
import sys
import scipy
import numpy
import scipy.special

# box dimensions (i.e. "extent" in ASPECT input)
xmin=0;xmax=11600.e3;
ymin=0;ymax=2900.e3;
# number of cells, (i.e. "number of repetitions" in ASPECT input)
xnum=5800
ynum=1450

y_crust = 7.5e3; 
slab_dip = 70.;

x_SP  = 4000.e3;
x_MP  = 2000.e3
x_OP  = 2000.e3 
x_gap = 0.5*(xmax - (x_SP + x_MP + x_OP))
gap_left = 0.75 * x_gap
gap_right = 1.25 * x_gap


depth_notch  = 250e3;
radius_outer = 250e3;
slab_dip = 70.;
OPthick = 100.e3;
weak_zone_width = 200e3

No_nodes= (xnum + 1) * (ynum + 1)
C=numpy.zeros([No_nodes,6],float)

SPcoredepth = 30.e3
SPcorethick = 25.e3

ind=0

for j in range(ynum + 1): 
	for i in range(xnum + 1):

		x = xmin + i * ((xmax - xmin)/xnum)
		y = ymin + j * ((ymax - ymin)/ynum) 

		C[ind,0] = x
		C[ind,1] = y

		# flat lying crust
		if x > (gap_left) and x <= (gap_left + x_SP - 0.5*radius_outer) and y > (ymax - y_crust):
			C[ind,2]=1	
		elif x > (gap_left + x_SP - 0.5*radius_outer + 100e3) and x < (gap_left + x_SP - 0.5*radius_outer + 500e3):
			tapered_crust = y_crust * (x - (gap_left + x_SP - 0.5*radius_outer + 100e3))/((gap_left + x_SP - 0.5*radius_outer + 500e3) - (gap_left + x_SP - 0.5*radius_outer + 100e3))
			if y > (ymax - tapered_crust): 
				C[ind,5]=1
		elif x >= (gap_left + x_SP - 0.5*radius_outer + 500e3) and x <= (gap_left + x_SP + x_MP) and y > (ymax - y_crust): 
				C[ind,5]=1
		elif x > (gap_left + x_SP + x_MP) and x < (xmax - gap_right) and y > (ymax - y_crust):
			C[ind,2]=1

		if x > (gap_left + x_SP - 0.5*radius_outer) and x <= (gap_left + x_SP + 0.5*radius_outer):
			x1 = gap_left + x_SP - 0.5*radius_outer; 
			y1 = ymax - radius_outer;
			if ((x-x1)**2 + (y-y1)**2) < radius_outer**2 and ((x-x1)**2 + (y-y1)**2) >= (radius_outer-y_crust)**2 and y > (ymax - depth_notch): 
				angle=numpy.arctan((y-y1)/(x-x1));
				if angle > numpy.radians(90. - slab_dip):
				   C[ind,2]=1

		if x > (gap_left + x_SP - 0.5*radius_outer) and x <= (gap_left + x_SP + 0.5*radius_outer):
			x1 = gap_left + x_SP - 0.5*radius_outer; 
			y1 = ymax - radius_outer;
			if ((x-x1)**2 + (y-y1)**2) >= radius_outer**2 and y > (ymax - OPthick) and y <= (ymax - y_crust): 
				C[ind,3]= 1

		if  x > (gap_left + x_SP + 0.5*radius_outer) and x < (gap_left + x_SP + x_MP) and y > (ymax - OPthick) and y <= (ymax - y_crust): 
				C[ind,3]= 1
		if x >= (gap_left + x_SP + x_MP) and x < (gap_left + x_SP + x_MP + weak_zone_width) and y > (ymax - OPthick) and y <= (ymax - y_crust): 
				C[ind,4]= 1

		ind=ind+1;

# write to file
f= open("text_files/compnotch_thick7.5km_SP4000kmMP2000kmOP2000km_BothPlatesCrustAndMP_DistinctCrustMPTapered.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(ynum+1)))
f.write("# Columns: x y composition1 composition2 composition3 composition4\n")
for k in range(0,ind):
	f.write("%.6f %.6f %.3f %.3f %.3f %.3f\n" % (C[k,0],C[k,1],C[k,2],C[k,3],C[k,4],C[k,5]))
f.close() 

