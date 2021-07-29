# Simulation of Trypsin

The necessary scripts and files are included to run GCMC/MD, GCNCMC/MD and water hopping simulatons on the equilibrated trypsin system (PDB code: 5MO2). Scripts are also included for the post-processing of these trajectories, as well as for analysis done both by clustering and electron density calculations.

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
