#!/bin/python3
import sys, os
import MDAnalysis as mda

iteration = 5
x, y, z = 11.5823, 11.5823, 17.5452 # in nanometers

os.system(f"insane -f CG-MD-{iteration}-solu-memb-initial.pdb -o CG-MD-{iteration}-solvated.gro -p topol.top -x {x} -y {y} -z {z} -sol W -salt 0.15 -center")
os.system("rm -v \#*")
