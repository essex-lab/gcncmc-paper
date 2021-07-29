# Analysis of simulations

The necessary scripts are included to perform two types of analyses:
1. Hierarchal based clustering analysis of the water locations 
2. Generation of electron density maps that can then be compared with those obtained experimentally

To perform a clustering analysis a co-ordinate file and post-processed trajectory are required:
```commandline
python cluster.py --coords <co-ordinate file> --traj <trajectory>
```

To calculate an electron density map the trajectory needs to be first aligned to the crystal structure and the first frame needs to be extracted. The number of frames in the trajectory needs to be entered as an argument. The Computational Crystallography Toolbox (CCTBX) is also required. This can be installed via conda via ```conda install -c conda-forge cctbx``` or via their github: https://github.com/cctbx/cctbx_project. Once installed, the following commands need to be run:
```commandline
python align_to_crystal_structure.py --coords <co-ordinate file> --traj <trajectory>
python get_first_frame.py --coords <co-ordinate file>
python xtraj.py traj=aligned_to_crystal.dcd top=first.pdb first=0 last=<number of frames>
```

The electron density calculation generates an ```fcalc.mtz``` file which can then be compared to the experimental electron density maps using software such as COOT or Phenix.
