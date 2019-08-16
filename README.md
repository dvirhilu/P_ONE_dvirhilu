# P_ONE_dvirhilu
Simulation and analysis of geometries and detector designs for Pacific Ocean Neutrino Explorer (P-ONE). P-ONE is a 200 optical module planned predecessor to a large scale neutrino telescope in the Cascadia Basin. 

Currently, The repository heavily relies on software made by the group for IceCube, a neutrino detector located in the South Pole. The repository does not contain said software, and to run the code in this repository, downloading and building IceCube software is essential. All simulation scripts were run on ComputeCanada's Cedar or University of Alberta's Illume computing clusters, and therefore all paths in simulation scripts refer to paths in the distributed system. All other scripts were run locally and have paths referencing local directories. It is also important to add P_ONE_dvirhilu/src to PYTHONPATH since some of the scripts reference modules from this directory.

Below is an outline of the file structure in the repository.

## src
This directory contains all source code used for the project.

### src/experimentModelCode
This directory contains any code used to make and analyze files describing the experiment model. This includes both the medium property files used by CLSim for photon propagation, and the DOM model files used for simulating hit detection. This directory contains the following files:

* FunctionClasses.py - A module containing classes representing different function types. At the moment only has classes representing polynomials (any degree) and functions from tables, as no other functions were needed.

* compareDOMModels.py - A script that takes DOM model files and produces plots detailing the DOM acceptance and angular acceptance
