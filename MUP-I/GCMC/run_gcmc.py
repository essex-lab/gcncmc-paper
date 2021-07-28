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

pdb = PDBFile('../MUP_crystal_equil.pdb')

# Add ghost waters
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology, pdb.positions,
                                                             n=15, pdb='MUP-prod-ghosts.pdb')

# Create system
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml', '../F09.xml')
system = forcefield.createSystem(pdb.topology, 
                                 nonbondedMethod=PME, 
                                 nonbondedCutoff=12*angstrom,
                                 switchDistance=10*angstrom, 
                                 constraints=HBonds)

# Add restraints to all heavy protein/ligand atoms
#restraint_force = CustomExternalForce('k*(((x-x0)^2)+((y-y0)^2)+((z-z0)^2))')
#restraint_force.addGlobalParameter('k', 10.0*kilocalories_per_mole/angstroms**2)
#restraint_force.addPerParticleParameter('x0')
#restraint_force.addPerParticleParameter('y0')
#restraint_force.addPerParticleParameter('z0')
#
#for residue in pdb.topology.residues():
#    if residue.name != 'HOH' and residue.name != 'Na+' and residue.name != 'Cl-':
#        for atom in residue.atoms():
#            if not atom.name.startswith('H'):
#                restraint_force.addParticle(atom.index, pdb.positions[atom.index].in_units_of(nanometers))
#
#system.addForce(restraint_force)

# Define reference atoms for the GCMC sphere
ref_atoms = [{'name': 'CA', 'resname': 'LYS', 'resid': '56', 'chain': 0},
             {'name': 'CA', 'resname': 'LEU', 'resid': '117', 'chain': 0}]

# Define integrator
integrator = BAOABIntegrator(298*kelvin, 1.0/picosecond, 0.002*picosecond)

# Create GCMC sampler object
gcmc_mover = grand.samplers.StandardGCMCSphereSampler(system=system, 
                                                      topology=pdb.topology, 
                                                      temperature=298*kelvin, 
                                                      referenceAtoms=ref_atoms, 
                                                      sphereRadius=6*angstroms, 
                                                      ghostFile="MUP-prod-ghosts.txt", 
                                                      log="MUP-prod.log", 
                                                      rst="MUP-prod.rst7", 
                                                      dcd="MUP-prod.dcd", 
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
#gcmc_mover.deleteWatersInGCMCSphere() # Uncomment to start from a dry sphere

# Run GCMD production (100k GCMC moves over 2.5 ns - 50 moves every 500 fs)
for i in range(5000):
    simulation.step(args.steps)
    gcmc_mover.move(simulation.context, args.moves)
    gcmc_mover.report(simulation)

