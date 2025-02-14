#!/bin/python3
from openmm.app import *
from openmm import *
from openmm.unit import *
'''
Hardware Setup
'''
# platform = Platform.getPlatformByName("CUDA")
# properties = {"DeviceIndex": "0", "Precision": "single"}
'''
Simulation System Setup
'''
iteration = 1
path = "." 
prmtop = AmberPrmtopFile(f"{path}/AA-MD-{iteration}/step5_input.parm7")
inpcrd = AmberInpcrdFile(f"{path}/AA-MD-{iteration}/step5_input.rst7")
# modeller = Modeller(pdb.topology, pdb.positions)
# modeller.addHydrogens(forcefield)
'''
Membrane Embedding and/or System Solvation
'''
# modeller.addMembrane(forcefield, lipidType="POPC", minimumPadding=1*nanometer)
# modeller.addSolvent(forcefield, model="tip3p",
                # padding=1*nanometer,
                # boxSize=Vec3(7.3, 7.3, 7.3)*nanometers,
                # numAdded=5000, # both water and ions
                # positiveIon="Na+", negativeIon="Cl-", 
                # ionicStrength=Quantity(value=0.15,unit=molar),
                # neutralize=True)
'''
Simulation Settings
'''
system = prmtop.createSystem(
    nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer, constraints=HBonds, 
    rigidWater=True, ewaldErrorTolerance=0.0005, hydrogenMass=1.0*amu,
)
barostat = MonteCarloMembraneBarostat(1*atmosphere, 0*bar*nanometer, 310*kelvin, MonteCarloMembraneBarostat.XYIsotropic,MonteCarloMembraneBarostat.ZFree, 10)
system.addForce(barostat)
integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 2*femtoseconds)
integrator.setConstraintTolerance(0.00001)
simulation = Simulation(prmtop.topology, system, integrator) # , platform, properties)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
'''
Energy Minimization
'''
simulation.minimizeEnergy(maxIterations=250000, tolerance=1.0)
positions = simulation.context.getState(
    getPositions=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions()).getPositions()
PDBFile.writeFile(simulation.topology, positions, open(f"{path}/AA-MD-{iteration}-initial.pdb", "w"))
'''
Equilibration
'''
simulation.context.setVelocitiesToTemperature(310*kelvin)
simulation.step(500000) # 1ns
'''
Production Simulation
'''
nb_steps = 37500000
simulation.reporters.append(XTCReporter(f"{path}/AA-MD-{iteration}/aa-md-{iteration}.xtc", 1000))
# simulation.reporters.append(DCDReporter(f"{path}/AA-MD-{iteration}/aa-md-{iteration}.dcd", 1000))
simulation.reporters.append(StateDataReporter(
    f"{path}/AA-MD-{iteration}/aa-md-{iteration}.log", 1000,
    totalSteps=nb_steps, step=True, speed=True, remainingTime=True,
    potentialEnergy=True, totalEnergy=True, density=True, temperature=True,
    volume=True, separator="\t"))
simulation.reporters.append(CheckpointReporter(f"{path}/AA-MD-{iteration}/aa-md-{iteration}.chk", 10000))
simulation.currentStep = 0
simulation.step(nb_steps) # 75ns
simulation.saveState(f"{path}/AA-MD-{iteration}/aa-md-{iteration}.rst")

positions = simulation.context.getState(
    getPositions=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions()).getPositions()
PDBFile.writeFile(simulation.topology, positions, open(f"{path}/AA-MD-{iteration}-final.pdb", "w"))
