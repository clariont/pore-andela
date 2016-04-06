# read me

############################################
make_coll_cyl.py
Makes a giant hexagonal 2D sheet, then "cuts out" a section to use for the cylinder.
Usage instructions: 
    enter (approximate) cyl_r, cyl_h, and d (lattice spacing). 
    the actual cyl_r and cyl_h will be given in units of d, the python script prints them out.

    Notes:
    - script is weird because cylinder radius and lattice spacing are not independent parameters.
    - Final cylinder has the lattice spacing you input, but the cylinder radius and height will be a little off. 

    polymer lengths are controlled by 'fixlen'.


############################################
lammps in.files
Equilibration procedure: first run ideal chains, then run LJ chains, then start translocation.
Translocation occurs with a constant force on the colloid and the polymer brush tail particles.




