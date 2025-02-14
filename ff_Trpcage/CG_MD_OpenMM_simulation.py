from openmm.unit import *
from openmm import *
from openmm.app import *
# from openmmtools.integrators import *
import martini_openmm as martini
'''
Hardware Setup
'''
# platform = Platform.getPlatformByName("CUDA")
# properties = {"DeviceIndex": "0", "Precision": "double"}

iteration = 9
path = "."
'''
Simulation Settings
'''
conf = GromacsGroFile(f"{path}/CG-MD-{iteration}/CG-MD-{iteration}-solvated.gro")
box_vectors = conf.getPeriodicBoxVectors()
top = martini.MartiniTopFile(f"{path}/CG-MD-{iteration}/system.top", periodicBoxVectors=box_vectors, epsilon_r=15)
system = top.create_system(nonbonded_cutoff=1.1*nanometer)
system.addForce(MonteCarloBarostat(1*bar, 310*kelvin))
integrator = LangevinIntegrator(310*kelvin, 10/picosecond, 20*femtoseconds)
# integrator = VelocityVerletIntegrator(20*femtoseconds)
# integrator.addConstrainPositions()
# integrator.setRandomNumberSeed(0)
simulation = Simulation(top.topology, system, integrator) # , platform, properties)
simulation.context.setPositions(conf.getPositions())
'''
Energy Minimization
'''
simulation.minimizeEnergy(maxIterations=50000,tolerance=1.0)
positions = simulation.context.getState(
    getPositions=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions()).getPositions()
PDBFile.writeFile(simulation.topology, positions, open(f"{path}/CG-MD-{iteration}-initial.pdb", "w"))
'''
Equilibration
'''
simulation.context.setVelocitiesToTemperature(310*kelvin)
simulation.step(50000) # 1ns
'''
Production Simulation
'''
nb_steps = 5000000
simulation.reporters.append(XTCReporter(f"{path}/CG-MD-{iteration}/cg-md-{iteration}.xtc", 1000))
# simulation.reporters.append(DCDReporter(f"{path}/CG-MD-{iteration}/cg-md-{iteration}.dcd", 1000))
simulation.reporters.append(StateDataReporter(
    f"{path}/CG-MD-{iteration}/cg-md-{iteration}.log", 1000,
    totalSteps=nb_steps, step=True, speed=True, remainingTime=True,
    potentialEnergy=True, totalEnergy=True, density=True, temperature=True,
    volume=True, separator="\t"))
simulation.reporters.append(CheckpointReporter(f"{path}/CG-MD-{iteration}/cg-md-{iteration}.chk", 10000))
simulation.currentStep = 0
simulation.step(nb_steps) # 15ns
simulation.saveState(f"{path}/CG-MD-{iteration}/cg-md-{iteration}.rst")

positions = simulation.context.getState(
    getPositions=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions()).getPositions()
PDBFile.writeFile(simulation.topology, positions, open(f"{path}/CG-MD-{iteration}-final.pdb", "w"))
