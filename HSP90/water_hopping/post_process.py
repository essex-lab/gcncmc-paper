import MDAnalysis as mda
import numpy as np
import grand
import mdtraj

# Recentre traj
trj = grand.utils.recentre_traj(trajectory='output/HSP90-example.nc', resname='PHE', name='CA', resid=124)

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

# Calculate waters in sphere at each frame
u = mda.Universe('../HSP90_equil.pdb', 'HSP90-prod-final.dcd')
water_oxygens = u.select_atoms('resname HOH and name O')

u_sphere = mda.Universe('HSP90-sphere.pdb')
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
