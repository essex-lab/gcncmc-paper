# Simulation of Major Urinary Protein I

The necessary scripts and files are included to run GCMC/MD, GCNCMC/MD and water hopping simulatons on the equilibrated MUP-I system (PDB code: 1ZNK). Scripts are also included for the post-processing of these trajectories in order for analysis to be performed. There are input co-ordinate files for the four different conformations of MUP-I, as described in the manuscript, with both equilibrated and dry binding sites (with the exception of the dry state, where the equilibrated and dry states are the same).

For the GCMC/MD and GCNCMC/MD simulations, to apply positional restraints on the heavy atoms of the protein and ligand, lines 31-43 need to be uncommented in the python script. To apply the same restraints for the water hopping simulations, lines 30-31 in ```run_water_hopping.py``` and lines 43-45 in ```example.yaml``` need to be uncommented.

To run and post-process a GCMC/MD simulation (5000 iterations of 20 GCMC moves and 4 ps MD):
```commandline
cd GCMC/
python run_gcmc.py
python post_process.py
```

To run and post-process a GCNCMC/MD simulation (4000 iterations of 1 NCMC move (npert=99, nprop=50, switching time=7 ps) and 5 ps MD):
```commandline
cd GCNCMC/
python run_gcncmc.py
python post_process.py
```

To run and post-process a water hopping simulation (4000 iterations of 1 NCMC move (npert=249, nprop=20, switching time=10 ps) and 5 ps MD):
```commandline
cd water_hopping/
python run_water_hopping.py
python post_process.py
```
