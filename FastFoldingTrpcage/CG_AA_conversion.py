#!/bin/python3
import mstool
import sys, os

simType = "protein" # "protein" or "NA"
iteration = 9
'''
Must Include CRYST1 Information in the PDB file. 
'''
mstool.Backmap(f"CG-MD-{iteration}-solu-final.pdb")
os.system("scp workdir/step4_final.pdb " + f"AA-MD-{iteration+1}-solu-initial.pdb")
os.system("mv -v workdir/ " + f"CG-MD-{iteration}/")
os.system("rm -v \#*")
