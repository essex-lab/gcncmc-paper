from WaterTranslation import WaterTranslationMove
from blues.moves import MoveEngine
from blues.simulation import *
import json
from blues.settings import Settings

opt = Settings(example.yaml).asDict()
structure = opt['Structure']
#restart = parmed.amber.Rst7(args.rst7)
#structure.positions = restart.positions
#structure.box = restart.box
print(json.dumps(opt, sort_keys=True, indent=2, skipkeys=True, default=str))


#Select move type
ligand = WaterTranslationMove(structure, water_name='WAT')

#Iniitialize object that selects movestep
ligand_mover = MoveEngine(ligand)

#Generate the openmm.Systems outside SimulationFactory to allow modifications
systems = SystemFactory(structure, ligand.atom_indices, opt['system'])

# Uncomment this section freeze atoms in the system
# Used for MD or BLUES production runs of the C60 buckyball and the water box with dividing graphene sheets
#systems.md = systems.freeze_atoms(structure, systems.md, **opt['freeze'])
#systems.alch = systems.freeze_atoms(structure, systems.alch, **opt['freeze'])

# Uncomment this section to restrain atoms in the MD and alchemical system
# Used for MD and BLUES production runs of MUP-1 and HSP90 protein-ligand systems to restrain all of the alpha carbons and the ligands
#systems.md = systems.restrain_positions(structure, systems.md, **opt['restraints'])
#systems.alch = systems.restrain_positions(structure, systems.alch, **opt['restraints'])

#Generate the OpenMM Simulations
simulations = SimulationFactory(systems, ligand_mover, opt['simulation'], opt['md_reporters'], opt['ncmc_reporters'])

# Energy minimize system
#simulations.md.minimizeEnergy(maxIterations=0) #minimize until convergence is reached
#state = simulations.md.context.getState(getPositions=True, getEnergy=True)
#print('Minimized energy = {}'.format(state.getPotentialEnergy().in_units_of(unit.kilocalorie_per_mole)))

# MD simulation - Uncomment for equilibration and MD production runs only
#simulations.md.step(opt['simulation']['nstepsMD'])

#BLUES Simulation - Only used for BLUES production runs
# Uncomment this section to run BLUES
blues = BLUESSimulation(simulations)
blues.run(**opt['simulation'])
