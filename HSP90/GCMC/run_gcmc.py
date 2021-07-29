from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
from sys import stdout
import numpy as np
import mdtraj
import argparse
import grand

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--moves', type=int, default=20, help='Number of GCMC moves per cycle')
parser.add_argument('-s', '--steps', type=int, default=2000, help='Number of MD steps per cycle')
args = parser.parse_args()

pdb = PDBFile('../input_files/HSP90_equil.pdb')

# Add ghost waters
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology, pdb.positions,
                                                             n=15, pdb='HSP90-prod-ghosts.pdb')

# Create system
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml', '../input_files/6G7.xml')
system = forcefield.createSystem(pdb.topology, 
                                 nonbondedMethod=PME, 
                                 nonbondedCutoff=12*angstrom,
                                 switchDistance=10*angstrom, 
                                 constraints=HBonds)

# Define reference atoms for the GCMC sphere
ref_atoms = [{'name': 'CA', 'resname': 'LEU', 'resid': '34', 'chain': 0},
             {'name': 'CA', 'resname': 'GLY', 'resid': '83', 'chain': 0}]

# Define integrator
integrator = BAOABIntegrator(298*kelvin, 1.0/picosecond, 0.002*picosecond)

# Create GCMC sampler object
gcmc_mover = grand.samplers.StandardGCMCSphereSampler(system=system, 
                                                      topology=pdb.topology, 
                                                      temperature=298*kelvin, 
                                                      referenceAtoms=ref_atoms, 
                                                      sphereRadius=6*angstroms, 
                                                      ghostFile="HSP90-prod-ghosts.txt", 
                                                      log="HSP90-prod.log", 
                                                      rst="HSP90-prod.rst7", 
                                                      dcd="HSP90-prod.dcd", 
                                                      overwrite=False)

# Define platform and precision
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Create simulation object and set positions, velocities and box vectors
simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(298*kelvin)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Initialise the GCMC mover
gcmc_mover.initialise(simulation.context, ghosts)

# Run GCMD production (100k GCMC moves over 2.5 ns - 50 moves every 500 fs)
for i in range(5000):
    simulation.step(args.steps)
    gcmc_mover.move(simulation.context, args.moves)
    gcmc_mover.report(simulation)

