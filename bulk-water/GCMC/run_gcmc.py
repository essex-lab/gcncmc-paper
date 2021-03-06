import time
import numpy as np
import argparse
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator

import grand

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--moves', type=int, help='Number of instantaneous GCMC moves to run per cycle')
parser.add_argument('-s', '--steps', type=int, help='Number of MD steps to run per cycle')
args = parser.parse_args()

# Load in a water box PDB...
pdb = PDBFile('../water_box-eq.pdb')

# Add ghosts
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology,
                                                             pdb.positions,
                                                             n=100,
                                                             pdb='water-ghosts.pdb')

ff = ForceField('tip3p.xml')

# Create a system using the pdb.topology generated by add_ghosts()
system = ff.createSystem(pdb.topology,
                         nonbondedMethod=PME,
                         nonbondedCutoff=12.0*angstroms,
                         switchDistance=10.0*angstroms,
                         constraints=HBonds)

# Langevin integrator
integrator = BAOABIntegrator(298*kelvin, 1.0/picosecond, 0.002*picoseconds)

# Create GCMC sampler object
gcmc_mover = grand.samplers.StandardGCMCSystemSampler(system=system,
                                                      topology=pdb.topology,
                                                      temperature=298*kelvin,
                                                      excessChemicalPotential=-6.09*kilocalories_per_mole,
                                                      standardVolume=30.345*angstroms**3,
                                                      boxVectors=np.array(pdb.topology.getPeriodicBoxVectors()),
                                                      log=f'{args.moves}moves-{args.steps}steps.log',
                                                      ghostFile=f"ghosts-{args.moves}moves-{args.steps}steps.txt",
                                                      overwrite=False)

# Define platform and precision
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Create Simulation object and set positions, velocities and box vectors
simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(298*kelvin)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Initialise GCMC mover
gcmc_mover.initialise(simulation.context, ghosts)

# Run for 12 hours and write out every move
start = time.time()

limit = 12 * 60 * 60

for i in range(10000000000000):
    simulation.step(args.steps)

    gcmc_mover.move(simulation.context, args.moves)

    gcmc_mover.report(simulation)

    if time.time() >= start + limit:
        break

