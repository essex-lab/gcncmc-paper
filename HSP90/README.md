# Simulation of Heat Shock Protein 90

The necessary scripts and files are included to run GCMC/MD, GCNCMC/MD and water hopping simulatons on the equilibrated HSP90 system (PDB code: 5J64). Scripts are also included for the post-processing of these trajectories, as well as for analysis done both by clustering and electron density calculations.

To run and post-process GCMC/MD simulation:
```commandline
cd GCMC/
python run_gcmc.py --moves 20 --steps 2000
python post_process.py
```
