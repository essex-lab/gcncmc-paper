import MDAnalysis as mda
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--coords', help='Input co-ordinate file (PDB)')
args = parser.parse_args()

u = mda.Universe(args.coords, 'aligned_to_crystal.dcd')

atoms = u.select_atoms('all')

with mda.Writer('first.pdb', atoms.n_atoms) as f:
    for ts in u.trajectory:
        if ts.frame == 0:
            f.write(atoms)
            break
