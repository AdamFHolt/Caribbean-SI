#!/usr/bin/env python 
#cp'd from make_tempSPandOP_halfspace80Ma25Ma_rad200km_BigBox_SmallerOP.py

import sys
import numpy as np
import scipy
import scipy.special

# box dimensions (i.e. "extent" in ASPECT input)
xmin=0;xmax=11600.e3;
ymin=0;ymax=2900.e3;
# number of cells, (i.e. "number of repetitions" in ASPECT input)
xnum=5800
ynum=1450

x_SP  = 4000.e3;
x_MP  = 2000.e3
x_OP  = 2000.e3 
x_gap = 0.5*(xmax - (x_SP + x_MP + x_OP))
gap_left = 0.75 * x_gap
gap_right = 1.25 * x_gap



depth_notch  = 250e3;
radius_outer = 250e3;
Tmax = 1573.; Tmin = 273.;
Tcutoff = 1325.;
slab_dip = 70.;

age_ma=60;
age=age_ma*1e6*365*24*60*60;
age_op_forearc_ma=20.
age_op_forearc=age_op_forearc_ma*1e6*365*24*60*60;
age_op_ma=100;
age_op=age_op_ma*1e6*365*24*60*60;
k = 1e-6

No_nodes= (xnum + 1) * (ynum + 1)
T=np.zeros([No_nodes,3],float)
 
ind=0

for j in range(ynum + 1): 

	for i in range(xnum + 1):

		x = xmin + i * ((xmax - xmin)/xnum)
		y = ymin + j * ((ymax - ymin)/ynum) 

		T[ind,0] = x
		T[ind,1] = y
		T[ind,2] = Tmax

		if x > (gap_left) and x <= (gap_left + x_SP - 0.5*radius_outer):
			ridge_dist = x_SP - 0.5*radius_outer
			age_ridge = (x - gap_left) * (age/ridge_dist)
			erf_term=(ymax-y)/(2*np.sqrt(k*age_ridge))
			T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))

		elif x > (gap_left + x_SP + 0.5*radius_outer) and x <= (gap_left + x_SP + x_MP):
			if x < (gap_left + x_SP + 500.e3):
				age_tmp = age_op_forearc + (age_op - age_op_forearc)*(x-(gap_left + x_SP - 0.5*radius_outer))/((gap_left + x_SP + 500.e3)-(gap_left + x_SP - 0.5*radius_outer))
				erf_term=(ymax-y)/(2*np.sqrt(k*age_tmp))
				T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))
			else:
				erf_term=(ymax-y)/(2*np.sqrt(k*age_op))
				T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))

		elif x > (gap_left + x_SP + x_MP) and x < (xmax - gap_right):
			ridge_dist = (xmax - gap_right) - (gap_left + x_SP + x_MP)
			age_ridge = ((xmax - gap_right)-x) * (age_op/ridge_dist)
			erf_term=(ymax-y)/(2*np.sqrt(k*age_ridge))
			T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))
				
		if x > (gap_left + x_SP - 0.5*radius_outer) and x <= (gap_left + x_SP + 0.5*radius_outer):
			x1 = gap_left + x_SP - 0.5*radius_outer; 
			y1 = ymax - radius_outer;
			if ((x-x1)**2 + (y-y1)**2) < radius_outer**2 and y > (ymax - depth_notch): 
				angle=np.arctan((y-y1)/(x-x1));
				if angle > np.radians(90. - slab_dip):
					ynotch = radius_outer - np.sqrt((x-x1)**2 + (y-y1)**2)
					erf_term=(ynotch)/(2*np.sqrt(k*age))
					T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))
			elif ((x-x1)**2 + (y-y1)**2) >= radius_outer**2 and y > (ymax - depth_notch): 
				age_tmp = age_op_forearc + (age_op - age_op_forearc)*(x-(gap_left + x_SP - 0.5*radius_outer))/((gap_left + x_SP + 500.e3)-(gap_left + x_SP - 0.5*radius_outer))
				erf_term=(ymax-y)/(2*np.sqrt(k*age_tmp))
				T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))

		if T[ind,2] > Tcutoff:
		   T[ind,2] = Tmax

		ind=ind+1;

# write to file
f= open("text_files/tempSPandOP_halfspace60Ma100Ma_SP4000kmMP2000kmOP2000km_rad200km_PlumeSubd_NoBaseBL_TaperedForearc.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(ynum+1)))
f.write("# Columns: x y temperature\n")
for k in range(0,ind):
	f.write("%.6f %.6f %.6f\n" % (T[k,0],T[k,1],T[k,2]))
f.close()

