# read me

############################################
**make_coll_cyl.py**
Makes a giant hexagonal 2D sheet, then "cuts out" a section to use for the cylinder.  The hexagonal points are
the grafted atoms.  Polymer beads are grown from the grafted atoms.


**Usage instructions:**


enter (approximate) cyl_r, cyl_h, and d (lattice spacing).  the final radius will be in units of 'd'.
if cyl_r = 5 and d = 2, the actual cylinder radius will be approximately 10. (explained in Notes).
the python script prints out the actual radius, height, and grafting density.

Notes:
- script is weird because cylinder radius and lattice spacing are not independent parameters.
- Final cylinder has the lattice spacing you input, but the cylinder radius and height will be a little off. 
- polymer lengths are controlled by 'fixlen'.


############################################
**lammps in.files**
Equilibration procedure: first run ideal chains, then run LJ chains, then start translocation.  each process has its own in.file.

Translocation occurs with a constant force on the colloid and the polymer brush tail particles.




