# Iterative Multiscale Molecular Dynamics (iMMD) in OpenMM: Accelerating Conformational Sampling of Biomolecular Systems by Iterating All-Atom and Coarse-Grained Molecular Dynamics Simulations
We combine the strengths of all-atom (AA) and coarse-grained (CG) MD simulations to develop the iterative multiscale MD (iMMD) simulation workflow in OpenMM by iterating between AA and CG MD simulations to enhance the sampling of biomolecular conformations. The AA-CG-AA iterations are repeated over multiple cycles, facilitating the accelerated sampling of biomolecular conformations at a CG level, while the finer atomistic interactions are refined with AA simulators.

The current implementation of iMMD is in the OpenMM simulation package but can be easily adopted for GROMACS. The initial input for iMMD can be a simple atomistic PDB structure of the protein of interest if the simulation system is mere solution or only involves common membrane phospholipid molecules, including DLPC, DLPE, DMPC, DOPC, DPPC, POPC, and POPE. iMMD in OpenMM can then automatically embed the protein in the phospholipid membrane of choice and/or solvate the simulation system with water and ions. If the simulation system involves non-phospholipid molecules, the user should set up simulation systems outside iMMD and provide the topology and coordinate files to start iMMD simulations.

The all-atom (AA) MD simulation is then carried out, using the CHARMM36m as the default force field parameter set. The final frame of the AA simulation is converted to its coarse-grained (CG) representation using the Martinize2 framework for protein and backward python script for the other components of the simulation system. The CG MD simulation is then performed to facilitate the transitions from one conformational state to another, with MARTINI3 as the default force field parameter set. The final frame of the CG simulation is converted back to its AA representation using the mstool python package for the protein, common lipids, water, and ions and backward python script for the uncommon lipid molecules. The AA-CG-AA cycles are repeated to continue the iMMD simulation. 
![Figure-1](https://github.com/user-attachments/assets/a6c95652-9458-4bf4-a13b-78a0c4bd2d3e)

The following *python* packages must be installed to perform iMMD simulations:
* OpenMM: https://openmm.org/
* Martinize2 and Vermouth: https://github.com/marrink-lab/vermouth-martinize
* DSSP: https://anaconda.org/salilab/dssp
* Insane: https://github.com/Tsjerk/Insane
* Martini in OpenMM: https://github.com/maccallumlab/martini_openmm
* mstool: https://mstool.readthedocs.io/en/0.1.0/installation.html
* backward: https://github.com/Tsjerk/MartiniTools/blob/master/backward.py and https://github.com/Tsjerk/MartiniTools/tree/master/Mapping

We demonstrate the enhanced sampling ability of iMMD over multiple AA-CG-AA iterations on four representative systems, including two globular proteins (fast-folding variant of Trpcage and Z-matrix protein of Mammarenavirus lassaense (LASV)) and two membrane proteins (Torpedo nicotinic acetylcholine receptor (nAChR) and REGN7663 Fab binding to the CXCR4 receptors in a heterogeneous cholesterol/phosphatidylcholine membrane lipid bilayer). In particular, iMMD captures the folding of the fast-folding Trpcage and Z-protein of LASV starting from their extended conformations as well as the binding of the REGN7663 Fab binding to the CXCR4s and dimerization of the CXCR4 receptors within the heterogenous membrane, within few AA-CG-AA iterations. Example folders containing the input files for performing iMMD simulations of these four simulation systems are included in this repository. It is recommended to run iMMD in OpenMM on GPUs to achieve the best simulation speeds.

# Funding
The workflow was developed under the supports by the National Institute of Allergy and Infectious Diseases of the National Institutes of Health and by the Duke Center for HIV Structural Biology, grant number U54-AI170752-01, and Sandia National Grand Challenges Funding, CAPSIID Project. This work was performed at the Los Alamos National Laboratory, which is operated by Triad National Security, LLC, for the National Nuclear Security Administration of the U.S. Department of Energy (contract 89233218CNA000001).  

# Reference
Do, H.N. and Gnanakaran, S. (2025) Iterative Multiscale Molecular Dynamics: Accelerating Conformational Sampling of Biomolecular Systems by Iterating All-Atom and Coarse-Grained Molecular Dynamics Simulations. BioRxiv. https://doi.org/10.1101/2025.02.16.638568 




