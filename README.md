# P_ONE_dvirhilu
Simulation and analysis of geometries and detector designs for Pacific Ocean Neutrino Explorer (P-ONE). P-ONE is a 200 optical module planned predecessor to a large scale neutrino telescope in the Cascadia Basin. 

Currently, The repository heavily relies on software made by the group for IceCube, a neutrino detector located in the South Pole. The repository does not contain said software, and to run the code in this repository, downloading and building IceCube software is essential. All simulation scripts were run on ComputeCanada's Cedar or University of Alberta's Illume computing clusters, and therefore all paths in simulation scripts refer to paths in the distributed system. All other scripts were run locally and have paths referencing local directories. It is also important to add P_ONE_dvirhilu/src to PYTHONPATH since some of the scripts reference modules from this directory.

Below is an outline of the file structure in the repository.

## src
This directory contains all source code used for the project.

### src/experimentModelCode
This directory contains any code used to make and analyze files describing the experiment model. This includes both the medium property files used by CLSim for photon propagation, and the DOM model files used for simulating hit detection. This directory contains the following files:

* FunctionClasses.py - A module containing classes representing different function types. At the moment only has classes representing polynomials (any degree) and functions from tables, as no other functions were needed.

* compareDOMModels.py - Takes DOM model files and produces plots detailing the DOM acceptance and angular acceptance

* makeDOMModelFiles.py - Creates the files detailing the DOM's acceptance vs. wavelength and the DOM's angular acceptance. These files are DOMEfficiency.dat and AngularAcceptance.dat and are contained in P_ONE_dvirhilu/DOMCharacteristics/.

* makeMedium.py - Creates the medium property files that are inputted into CLSim. These files include cfg.txt, icemodel.dat, icemodel.par (named icemodel because that's the name CLSim takes, but describe water properties), and are located in P_ONE_dvirhilu/propagationMediumModels

### src/gcdScripts
This directory contains all the scripts used to make the different geometries (GCD I3 files). The only geometry files not produced by scripts in this directory are partialDenseGeo, and are produced by P_ONE_dvirhilu/src/SimAnalysis/frameNumGeoFilter.py . GCD files are too big to be placed in a github repository, but are located in the desktop that was used during the workterm. All gcd files are in /home/dvir/workFolder/I3Files/gcd/.Sub-directories in the gcd directory are all named according to the type of geometry they represent, as aligned with the scripts detailed in src/gcdScripts. To get an understanding of what the different geometries look like, use the steamshovel event viewer. To understand what the names of the GCD files mean, visit the script that produced them. All gcd scripts have input paramaters that can change the specific details of that geometry type.

* gcdHelpers.py - contains helper functions used across many of the GCD producing scripts

* generateCubeGeometry.py - creates a cube shaped geometry

* generateDenseGeometry.py - creates the large, dense geometry used for geometry optimization

* generateHorizontalGeometryP_ONE.py - generates the Horizontal geometry proposed for P-ONE

* getCDFramesData.py - takes the much of the calibration and detector status data from IceCube and averages it to input it into the gcd files

* makeTestString.py - creates a single, horizontal string



