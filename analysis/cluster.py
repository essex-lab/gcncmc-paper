import grand
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--coords', help='Input co-ordinate file (PDB)')
parser.add_argument('-t', '--traj', help='Post-processed trajectory')
args = parser.parse_args()

if 'MUP' in args.coords:
    ref_atoms = [{'name': 'CA', 'resname': 'LYS', 'resid': '56', 'chain': 0},
                 {'name': 'CA', 'resname': 'LEU', 'resid': '117', 'chain': 0}]
elif 'HSP90' in args.coords:
    ref_atoms = [{'name': 'CA', 'resname': 'LEU', 'resid': '34', 'chain': 0},
                 {'name': 'CA', 'resname': 'GLY', 'resid': '83', 'chain': 0}]
elif 'trypsin' in args.coords:
    ref_atoms = [{'name': 'CA', 'resname': 'GLY', 'resid': '204', 'chain': 0},
                 {'name': 'CA', 'resname': 'ALA', 'resid': '198', 'chain': 0}]
else:
    print('Unable to determine reference atoms')

grand.utils.cluster_waters(topology=args.coords, trajectory=args.traj, sphere_radius=6.0, ref_atoms=ref_atoms)
