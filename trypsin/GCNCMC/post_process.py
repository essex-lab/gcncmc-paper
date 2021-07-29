import grand
import mdtraj

# Shift ghost waters
trj = grand.utils.shift_ghost_waters(ghost_file='trypsin-prod-ghosts.txt',
                                     topology='trypsin-prod-ghosts.pdb',
                                     trajectory='trypsin-prod.dcd')

# Recentre traj
trj = grand.utils.recentre_traj(t=trj, resname='GLY', name='CA', resid=189)

# Align traj
grand.utils.align_traj(t=trj, output='trypsin-prod-final.dcd')

# Write sphere
ref_atoms = [{'name': 'CA', 'resname': 'GLY', 'resid': '204', 'chain': 0},
             {'name': 'CA', 'resname': 'ALA', 'resid': '198', 'chain': 0}]

grand.utils.write_sphere_traj(radius=6,
                              ref_atoms=ref_atoms,
                              topology='trypsin-prod-ghosts.pdb',
                              trajectory='trypsin-prod-final.dcd',
                              output='trypsin-sphere.pdb',
                              initial_frame=True)
