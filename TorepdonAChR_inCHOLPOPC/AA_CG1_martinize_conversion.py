#!/bin/python3
import sys, os
import MDAnalysis as mda

iteration = 5
dssp = "/home/hungd238/miniconda3/envs/multiscaleEnv/bin/mkdssp"

os.system(f"martinize2 -f AA-MD-{iteration}-solu-final.pdb -x CG-MD-{iteration}-solu-initial.pdb -o system.top -resid input -p backbone -pf 2500 -dssp {dssp} -from charmm -ff martini3001 -elastic -ef 1500.0 -el 0 -eu 1.2 -ea 0 -ep 1.0 -eunit chain -scfix -cys auto -nt -maxwarn 1")

os.system(f"python backward.py -f AA-MD-{iteration}-POPC-final.pdb -o CG-MD-{iteration}-POPC-initial.gro -from charmm36 -to martini")
univ = mda.Universe(f"CG-MD-{iteration}-POPC-initial.gro")
with mda.Writer(f"CG-MD-{iteration}-POPC-initial.pdb") as pdb:
    pdb.write(univ)

os.system(f"python backward.py -f AA-MD-{iteration}-CHOL-final.pdb -o CG-MD-{iteration}-CHOL-initial.gro -from amber -to martini")
univ = mda.Universe(f"CG-MD-{iteration}-CHOL-initial.gro")
with mda.Writer(f"CG-MD-{iteration}-CHOL-initial.pdb") as pdb:
    pdb.write(univ)
os.system("rm -v \#*")
