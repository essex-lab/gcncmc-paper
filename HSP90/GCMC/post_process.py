import grand
import MDAnalysis as mda
from MDAnalysis.analysis import align
import mdtraj

# Shift ghost waters
trj = grand.utils.shift_ghost_waters(ghost_file='HSP90-prod-ghosts.txt',
                                     topology='HSP90-prod-ghosts.pdb',
                                     trajectory='HSP90-prod.dcd')

# Recentre traj
trj = grand.utils.recentre_traj(t=trj, resname='PHE', name='CA', resid=124)

# Align traj
grand.utils.align_traj(t=trj, output='HSP90-prod-final.dcd')

# Write sphere
ref_atoms = [{'name': 'CA', 'resname': 'LEU', 'resid': '34', 'chain': 0},
             {'name': 'CA', 'resname': 'GLY', 'resid': '83', 'chain': 0}]

grand.utils.write_sphere_traj(radius=6,
                              ref_atoms=ref_atoms,
                              topology='HSP90-prod-ghosts.pdb',
                              trajectory='HSP90-prod-final.dcd',
                              output='HSP90-sphere.pdb',
                              initial_frame=True)
