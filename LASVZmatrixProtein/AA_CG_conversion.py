#!/bin/python3
import sys, os

iteration = 6
x, y, z = 8.7671, 8.7671, 8.7671 # in nanometers
dssp = "/home/hungd238/miniconda3/envs/multiscaleEnv/bin/mkdssp"

os.system(f"martinize2 -f AA-MD-{iteration}-solu-final.pdb -x CG-MD-{iteration}-solu-initial.pdb -o system.top -resid input -p backbone -pf 2500 -dssp {dssp} -from charmm -ff martini3001 -elastic -ef 1500.0 -el 0 -eu 0.9 -ea 0 -ep 1.0 -eunit chain -scfix -cys auto -nt -maxwarn 1")

os.system(f"insane -f CG-MD-{iteration}-solu-initial.pdb -o CG-MD-{iteration}-solvated.gro -p system.top -x {x} -y {y} -z {z} -sol W -salt 0.15 -center") 
# -d 5.0 -pbc cubic -sol W -salt 0.15 -center")
os.system("rm -v \#*")
