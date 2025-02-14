#!/bin/python3
import MDAnalysis as mda
import mstool
import sys, os

simType = "protein" # "protein" or "NA"
iteration = 5
'''
Must Include CRYST1 Information in the PDB file. 
'''
mstool.Backmap(f"CG-MD-{iteration}-solu-final.pdb")
os.system("scp workdir/step4_final.pdb " + f"AA-MD-{iteration+1}-solu-initial.pdb")
os.system("mv -v workdir/ " + f"CG-MD-{iteration}/workdir-solu")

mstool.Backmap(AA=f"AA-MD-{iteration}-solu-final.pdb", structure=f"CG-MD-{iteration}-memb-final.pdb")
os.system("scp workdir/step4_final.pdb " + f"AA-MD-{iteration+1}-solu-memb-initial.pdb")
os.system("mv -v workdir/ " + f"CG-MD-{iteration}/")
os.system("rm -v \#*")
