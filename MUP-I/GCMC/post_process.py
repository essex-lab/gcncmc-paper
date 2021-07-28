import grand
import mdtraj

# Shift ghost waters
trj = grand.utils.shift_ghost_waters(ghost_file='MUP-prod-ghosts.txt',
                                     topology='MUP-prod-ghosts.pdb',
                                     trajectory='MUP-prod.dcd')

# Recentre traj
trj = grand.utils.recentre_traj(t=trj, resname='F09', name='C5', resid=1)

# Align traj
grand.utils.align_traj(t=trj, output='MUP-prod-final.dcd')

# Write sphere
ref_atoms = [{'name': 'CA', 'resname': 'LYS', 'resid': '56', 'chain': 0},
             {'name': 'CA', 'resname': 'LEU', 'resid': '117', 'chain': 0}]

grand.utils.write_sphere_traj(radius=6,
                              ref_atoms=ref_atoms,
                              topology='MUP-prod-ghosts.pdb',
                              trajectory='MUP-prod-final.dcd',
                              output='MUP-sphere.pdb',
                              initial_frame=True)
