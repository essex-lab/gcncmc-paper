# Simulation of Heat Shock Protein 90

The necessary scripts and files are included to run GCMC/MD, GCNCMC/MD and water hopping simulatons on the equilibrated HSP90 system (PDB code: 5J64). Scripts are also included for the post-processing of these trajectories, as well as for analysis done both by clustering and electron density calculations.

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
