#!/bin/python3
import sys, os

iteration = 9
x, y, z = 7.0789, 7.0789, 7.0789 # in nanometers
dssp = "/home/hungd238/miniconda3/envs/multiscaleEnv/bin/mkdssp"

os.system(f"martinize2 -f AA-MD-{iteration}-solu-final.pdb -x CG-MD-{iteration}-solu-initial.pdb -o system.top -resid input -p backbone -pf 2500 -dssp {dssp} -from charmm -ff martini3001 -elastic -ef 1500.0 -el 0 -eu 1.2 -ea 0 -ep 1.0 -eunit chain -scfix -cys auto -nt -maxwarn 1")

os.system("insane -f " + f"CG-MD-{iteration}-solu-initial.pdb -o CG-MD-{iteration}-solvated.gro" + f" -p system.top -x {x} -y {y} -z {z} -sol W -salt 0.15 -center")
os.system("rm -v \#*")
