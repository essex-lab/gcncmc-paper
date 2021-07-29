import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import align

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--coords', help='Input co-ordinate file (PDB)')
parser.add_argument('-t', '--traj', help='Post-process trajectory')
args = parser.parse_args()

u = mda.Universe(args.coords, args.traj)

u_ref = mda.Universe('../HSP90/input_files/HSP90_aligned_to_crystal.pdb')

alignment = align.AlignTraj(u, u_ref, filename="aligned_to_crystal.dcd", select="protein and name CA")
alignment.run()
