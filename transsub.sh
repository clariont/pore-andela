#!/bin/bash

lmpdir='/home/clarion/Programs/lmp_protein/src/lmp_serial'

python readprotmakecyl2.py


# Write new cylinder radius to LAMMPS input files:
myrad=$( cat myrad.dat )
head -n 6 in.ideal > temp
echo "variable cylr equal $myrad" >> temp
mylines=$( wc -l < in.ideal )
tail -n $(( mylines - 7 )) in.ideal >> temp
mv temp in.ideal

head -n 6 in.real1 > temp
echo "variable cylr equal $myrad" >> temp
mylines=$( wc -l < in.real1 )
tail -n $(( mylines - 7 )) in.real1 >> temp
mv temp in.real1

head -n 6 in.flow > temp
echo "variable cylr equal $myrad" >> temp
mylines=$( wc -l < in.flow )
tail -n $(( mylines - 7 )) in.flow >> temp
mv temp in.flow

head -n 6 in.trans > temp
echo "variable cylr equal $myrad" >> temp
mylines=$( wc -l < in.trans )
tail -n $(( mylines - 7 )) in.trans >> temp
mv temp in.trans


# Start trans:
time $lmpdir < in.ideal
cp ideal1.dat ideal.dat
time $lmpdir < in.real
cp real1.dat real.dat
time $lmpdir < in.flow
cp flow1.dat flow.dat
time $lmpdir < in.trans



