import MDAnalysis as mda
import numpy as np
import grand
import mdtraj

# Recentre traj
trj = grand.utils.recentre_traj(trajectory='output/MUP-example.nc', resname='F09', name='C5', resid=1)

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

# Calculate waters in sphere at each frame
u = mda.Universe('../MUP_crystal_equil.pdb', 'MUP-prod-final.dcd')
water_oxygens = u.select_atoms('resname HOH and name O')

u_sphere = mda.Universe('MUP-sphere.pdb')
sphere = u_sphere.select_atoms('resname SPH')

waters = []

for ts in range(len(u.trajectory)):
    u.trajectory[ts]
    u_sphere.trajectory[ts]

    x = 0

    for o in water_oxygens:
        if np.linalg.norm(o.position - sphere.positions[0]) < 6.0:
            x += 1

    waters.append(x)

with open('waters.txt', 'w') as f:
    for water in waters:
        f.write(f'{water}\n')
