#!/usr/bin/python

############################################################################################################
##############################
##  Date: 5 Apr 2016
##  Description: python script to make a pore
##		 
##		 
##
##  Changelog:
##  
##
##############################
############################################################################################################

import math
import random
import sys


#System Parameters:
nhex = 100					# don't need to change this!
cyl_r = 2.0					# Cylinder Radius (approximate)
cyl_h = 10					# Cylinder Height
xlo = -100					# box bounds
xhi = 100
ylo = -100
yhi = 100
zlo = -500
zhi = 500
d = 4					    	# lattice spacing

# Chain params:
fixlen = 20		# number of beads in a chain	
sigma = 1.0		# bead sigma 
Rg = sigma*0.5
bond_r = 1.2		# bond length for polymer beads


# Colloid:
ncoll = 1





##########################################################################################################
## Main
##########################################################################################################

cyl_r = cyl_r*d
cyl_h = cyl_h*d

cyl_circ = math.pi*2*cyl_r
circ1 = cyl_circ-d

# Make Hexagonal lattice, map it to a cylinder
npts = nhex*nhex*nhex - (nhex-1)*(nhex-1)*(nhex-1)
points = []
myShift = math.sqrt(3)*d/4.0  # we want an even number of circles for the cylinder
for i in range(nhex):
    y = math.sqrt(3)*i*d/2.0
    for j in range(2*nhex-1-i):
	x = (-(2*nhex-i-2)*d)/2.0+j*d;
	points.append([x,y+myShift,0])
	if (y != 0):
	    points.append([x,(-y+myShift),0])
realpts = []
xlowest = 0
xgreatest = 0
ylowest = 0
ygreatest = 0
zlowest = 0
zgreatest = 0
for pt in points:
    if (pt[0] > -circ1/2.0 and pt[0] < (circ1/2.0+d/2.0+0.000001) and pt[1] < cyl_h/2.0 and pt[1] > -cyl_h/2.0):
	if (pt[0] < xlowest):
	    xlowest = pt[0]
	if (pt[0] > xgreatest):
	    xgreatest = pt[0]
	if (pt[0] < ylowest):
	    ylowest = pt[1]
	if (pt[1] > ygreatest):
	    ygreatest = pt[1]
	realpts.append(pt)

real_r = xgreatest - xlowest + d/2.0
#print "xlowest, xgreatest:", xlowest, xgreatest
#print "ylowest, ygreatest:", ylowest, ygreatest
#print "cyl top, cyl bottom: ", (ylowest-myShift), (ygreatest+myShift)
real_h = ygreatest - ylowest
cyl_circ = real_r
rr = cyl_circ/2/math.pi


cyl_volume = real_h*rr
rho = len(realpts)/cyl_volume
print "cyl_radius:", rr
print "cyl_radius + 0.5*sigma (for wall):", rr+1*Rg
print "cyl h: ", real_h
print "grafting density: ", len(realpts)/(2*math.pi*rr*real_h)

# Put on circle:
points = []
cylpoints = []
for pt in realpts:
    x = rr*math.cos(pt[0]*(2*math.pi)/cyl_circ)
    y = rr*math.sin(pt[0]*(2*math.pi)/cyl_circ)
    z = pt[1]
    points.append(pt)
    cylpoints.append([x,y,z])
#    print x,y,z

nlig = len(cylpoints)
#print "nlig: ", nlig

liglens = []
ligadd = 0
for i in range(len(cylpoints)):
    liglen = fixlen
    liglens.append(int(liglen))
    ligadd = ligadd + int(liglen)


# Write LAMMPS configuration file
atoms = ligadd + ncoll
bonds = ligadd-nlig 
angles = (fixlen-2)*(nlig)
dihedrals = 0
impropers = 0
atomtypes = (2+1)
bondtypes = 1
angletypes = 1
dihedraltypes = 0
impropertypes = 0

g=open("cylligs.dat", "w")
g.write("LAMMPS Description\n\n")
g.write("\t"+str(atoms)+"\tatoms\n")
g.write("\t"+str(bonds)+"\tbonds\n")
g.write("\t"+str(angles)+"\tangles\n")
g.write("\t"+str(dihedrals)+"\tdihedrals\n")
g.write("\t"+str(impropers)+"\timpropers\n\n")
g.write("\t"+str(atomtypes)+"\tatom types\n")
g.write("\t"+str(bondtypes)+"\tbond types\n")
g.write("\t"+str(angletypes)+"\tangle types\n")
g.write("\t"+str(dihedraltypes)+"\tdihedral types\n")
g.write("\t"+str(impropertypes)+"\timproper types\n\n\n")

#Write Box Lengths:
g.write(str(xlo)+" "+str(xhi)+" xlo xhi\n")
g.write(str(ylo)+" "+str(yhi)+" ylo yhi\n")
g.write(str(zlo)+" "+str(zhi)+" zlo zhi\n\n\n")

#Write Atomic Masses:
l = 1.0
g.write("Masses\n\n")
g.write("\t1 "+str(l)+"\n")
g.write("\t2 "+str(l)+"\n")
g.write("\t3 "+str(l)+"\n")

#Write Atoms:
g.write("\nAtoms\n\n")

num = 1   # counts the total number of atoms
lig1pos = []
lig2pos = []
allligs = []
overlap = 0



# Write Chain Monomers
for i in range(nlig):
    g.write("\t"+str(num)+" "+str(num)+" 1 0 ")
    # Grafted atom
    x = cylpoints[i][0]
    y = cylpoints[i][1]
    z = cylpoints[i][2]
    g.write(str(x)+" "+str(y)+" "+str(z)+" 0 0 0\n")
    num = num + 1
    xg = x
    yg = y
    # Tail atoms
    # Zigzag monomers
    mag = math.sqrt(xg*xg+yg*yg)
    px=yg/mag			    # perpendicular to (xg,yg) and z-axis
    py=-xg/mag
    side=1
    for j in range(liglens[i]-1):
	myRandLen = random.randint(10,20)
	myRandLen = 15
	if (j < myRandLen):
	    x1 = xg-0.5*(j+1)*xg*bond_r/mag
	    y1 = yg-0.5*(j+1)*yg*bond_r/mag
	    if (z > 0):
		z1 = z - (j+1)*0.05*bond_r
	    else:
		z1 = z + (j+1)*0.05*bond_r
	    if (j%2 == 0):
		x1 = x1 - math.sqrt(3)/2.0*bond_r*px*side
		y1 = y1 - math.sqrt(3)/2.0*bond_r*py*side
		side = side*(-1)
	    else:
		x1 = x1+0.5*xg*bond_r/mag
		y1 = y1+0.5*yg*bond_r/mag
	    x0 = x1
	    y0 = y1
	    z0 = z1
	else:
	    # Random walk
	    evec = 0
	    while (evec == 0):
		ex = random.gauss(0,1)
		ey = random.gauss(0,1)
#		ez = random.gauss(0,1)
		ez = 0
		emag = math.sqrt(ex*ex+ey*ey+ez*ez)
		ex = ex/emag*bond_r*1
		ey = ey/emag*bond_r*1
		ez = ez/emag*bond_r*1
		x1 = x0 + ex
		y1 = y0 + ey
		z1 = z0 + ez
		if (((x1)*(x1)+(y1)*(y1) < rr*rr) and (z1*z1 < (cyl_h-2*Rg)*(cyl_h-2*Rg))):
		    evec = 1
		else:
		    evec = 0
	    x0 = x1
	    y0 = y1
	    z0 = z1

	g.write("\t"+str(num)+" "+str(num)+" 2 0 ")
	g.write(str(x1)+" "+str(y1)+" "+str(z1)+" 0 0 0\n")
	num = num + 1

# Write Colloid
# First translate the folded protein:
xcoll = 0
ycoll = 0
zcoll = -100
for i in range(ncoll):
    g.write("\t"+str(num)+" "+str(num)+" 3 0 ")
    g.write(str(xcoll)+" "+str(ycoll)+" "+str(zcoll)+" 0 0 0\n")
    zcoll = zcoll + 2
    num = num + 1
#    ctr = ctr + 1
 


# Write Bonds:
g.write("\nBonds\n\n")
countlig = 0
bondnum = 1
# Brushes
for i in range(nlig):
    g.write("\t"+str(bondnum)+" 1 "+str(1+i*(fixlen))+" "+str(2+i*(fixlen))+"\n")
    bondnum = bondnum + 1 
    # Bind head to first tail atom:
    for j in range(liglens[i]-2):
	g.write("\t"+str(bondnum)+" 1 "+str(i*(fixlen)+2+j)+" "+str(i*(fixlen)+3+j)+"\n")
	bondnum = bondnum + 1


# Write Angles:
g.write("\nAngles\n\n")
anglenum = 1
for i in range(nlig):
    for j in range(fixlen-2):
	g.write("\t"+str(anglenum)+" 1 "+str(i*(fixlen)+j+1)+" "+str(i*fixlen+j+2)+" "+str(i*fixlen+j+3)+"\n")
	anglenum = anglenum+1
    


# Print radius:
g=open("myrad.dat",'w')
g.write(str(cyl_r)+"\n")
g.close()

# Tcl script for VMD
g=open("setrad.tcl",'w')
g.write('set sel [atomselect top "all"]'+"\n")
g.write('$sel set radius'+" "+str(Rg)+"\n")
g.write('$sel delete'+"\n")
