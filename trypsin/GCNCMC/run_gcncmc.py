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
parser.add_argument('--pert', type=int, default=99, help='Number of perturbation steps per NCMC move')
parser.add_argument('--prop', type=int, default=50, help='Number of propagation steps between each perturbation')
args = parser.parse_args()

pdb = PDBFile('../input_files/trypsin_equil.pdb')

# Add ghost waters
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology, pdb.positions,
                                                             n=15, pdb='trypsin-prod-ghosts.pdb')

# Create system
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml', '../input_files/BEN.xml')
system = forcefield.createSystem(pdb.topology, 
                                 nonbondedMethod=PME, 
                                 nonbondedCutoff=12*angstrom,
                                 switchDistance=10*angstrom, 
                                 constraints=HBonds)

# Define reference atoms for the GCMC sphere
ref_atoms = [{'name': 'CA', 'resname': 'GLY', 'resid': '204', 'chain': 0},
             {'name': 'CA', 'resname': 'ALA', 'resid': '198', 'chain': 0}]

# Define integrator
integrator = BAOABIntegrator(298*kelvin, 1.0/picosecond, 0.002*picosecond)

# Create GCMC sampler object
gcncmc_mover = grand.samplers.NonequilibriumGCMCSphereSampler(system=system, 
                                                              topology=pdb.topology, 
                                                              temperature=298*kelvin, 
                                                              integrator=integrator,
                                                              nPertSteps=args.pert,
                                                              nPropStepsPerPert=args.prop,
                                                              referenceAtoms=ref_atoms, 
                                                              sphereRadius=6*angstroms, 
                                                              ghostFile="trypsin-prod-ghosts.txt", 
                                                              log="trypsin-prod.log", 
                                                              rst="trypsin-prod.rst7", 
                                                              dcd="trypsin-prod.dcd", 
                                                              overwrite=False)

# Define platform and precision
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Create simulation object and set positions, velocities and box vectors
simulation = Simulation(pdb.topology, system, gcncmc_mover.compound_integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(298*kelvin)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Initialise the GCMC mover
gcncmc_mover.initialise(simulation.context, ghosts)

# Run GCNCMC/MD production - 4000 iterations of 1 NCMC move and 5 ps MD
for i in range(4000):
    simulation.step(2500)
    gcncmc_mover.move(simulation.context, 1)
    gcncmc_mover.report(simulation)

