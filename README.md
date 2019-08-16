# P_ONE_dvirhilu
Simulation and analysis of geometries and detector designs for Pacific Ocean Neutrino Explorer (P-ONE). P-ONE is a 200 optical module planned predecessor to a large scale neutrino telescope in the Cascadia Basin. 

Currently, The repository heavily relies on software made by the group for IceCube, a neutrino detector located in the South Pole. The repository does not contain said software, and to run the code in this repository, downloading and building IceCube software is essential. All simulation scripts were run on ComputeCanada's Cedar or University of Alberta's Illume computing clusters, and therefore all paths in simulation scripts refer to paths in the distributed system. All other scripts were run locally and have paths referencing local directories. It is also important to add P_ONE_dvirhilu/src to PYTHONPATH since some of the scripts reference modules from this directory.

Below is an outline of the file structure in the repository.

## src
This directory contains all source code used for the project.

### src/experimentModelCode
This directory contains any code used to make and analyze files describing the experiment model. This includes both the medium property files used by CLSim for photon propagation, and the DOM model files used for simulating hit detection. This directory contains the following files:

* FunctionClasses.py - a module containing classes representing different function types. At the moment only has classes representing polynomials (any degree) and functions from tables, as no other functions were needed.

* compareDOMModels.py - takes DOM model files and produces plots detailing the DOM acceptance and angular acceptance

* makeDOMModelFiles.py - creates the files detailing the DOM's acceptance vs. wavelength and the DOM's angular acceptance. These files are DOMEfficiency.dat and AngularAcceptance.dat and are contained in P_ONE_dvirhilu/DOMCharacteristics/.

* makeMedium.py - creates the medium property files that are inputted into CLSim. These files include cfg.txt, icemodel.dat, icemodel.par (named icemodel because that's the name CLSim takes, but describe water properties), and are located in P_ONE_dvirhilu/propagationMediumModels

### src/gcdScripts
This directory contains all the scripts used to make the different geometries (GCD I3 files). The only geometry files not produced by scripts in this directory are partialDenseGeo, and are produced by P_ONE_dvirhilu/src/SimAnalysis/frameNumGeoFilter.py . GCD files are too big to be placed in a github repository, but are located in the desktop that was used during the workterm. All gcd files are in /home/dvir/workFolder/I3Files/gcd/.Sub-directories in the gcd directory are all named according to the type of geometry they represent, as aligned with the scripts detailed in src/gcdScripts. To get an understanding of what the different geometries look like, use the steamshovel event viewer. To understand what the names of the GCD files mean, visit the script that produced them. All gcd scripts have input paramaters that can change the specific details of that geometry type.

* gcdHelpers.py - contains helper functions used across many of the GCD producing scripts

* generateCubeGeometry.py - creates a cube shaped geometry

* generateDenseGeometry.py - creates the large, dense geometry used for geometry optimization

* generateHorizontalGeometryP_ONE.py - generates the Horizontal geometry proposed for P-ONE

* getCDFramesData.py - takes the much of the calibration and detector status data from IceCube and averages it to input it into the gcd files

* makeTestString.py - creates a single, horizontal string

### src/simAnalysis

This directory contains all the scripts used to analyze simulation results. 

* SimAnalysis.py - contains helper functions used across many of the simulation analysis scripts

* analyzeRecoQuality.py - Given a reconstruction I3File, the script analyzes how close the reconstructed particle is to the original track for each frame

* compareAeff.py - calculates the effective areas of I3Files outputted from simulations. It is important to use the I3File outputted directly from hit detection in order to not wrongfully count the number of events. Multiple I3Files can be inputted and compared together on the same plot.

* findNumHits.py - an early script used to compare the results of the IceCube hit detection simulation and the custom made script in order to test the custom simulation.

* frameNumGeoFilter.py - Takes the dense geometry and given information about the specific DOMs to pick out, calculates how many frames were retained and compares it against the initial tested geometry. If the frame number is comparable, it saves the geometry in /home/dvir/workFolder/I3Files/gcd/partialDenseGeo.

* getPhaseFunctionPlot.py - A simple script to compare different phase functions for photon scattering

* linefitReco.py - Takes simulation files after hit detection, filters using event selection criteria, and reconstructs the particle using the LineFit algorithm. It then saves this as a new I3 file in /home/dvir/workFolder/I3Files/linefitReco

* improvedTrackReco.py - Takes simulation files after hit detection, filters using event selection criteria, and reconstructs the particle using the hit time chi squared algorithm. It then saves this as a new I3 file in /home/dvir/workFolder/I3Files/linefitReco. Since this algorithm uses a minmizer, it outputs the minimizer results to MinimizerOutputs.txt 

* muonFluxDistributions.py - an early script that was used to analyze the distribution of mouns from a MuonGun simulation across different parameters

* nugenNeutrinoDistribution.py - an early script that was used to analyze the distribution of neutrinos from a Neutrino-Generator simulation across different parameters

* plotDOMTimeResiduals.py - plots the time residuals from a simulation file after photon propagation and compares it to data taken by STRAW. 

### src/simCode

Contains all scripts to run a full simulation with either Genie, MuonGun, or Neutrino-Generator. These files are outputted to /home/dvir/workFolder/I3Files and are sorted in folders by <generator-tool>/<simulation-step>/<gcd-type>/ . Step 1 refers to a file after particle generation, interaction, and propagation. Step 2 refers to a file after photon propagation. Step 3 refers to a file after hit detection.
  
NOTE: only the Neutrino-Generator scripts have the closest approach filter

* All simulation steps have to be submitted on a cluster

* On cedar, the submit file is labeled <name>_submit.sh . On illume, it is named <name>_condor.submit .

* On the cluster, only the submit file needs to be called. The submit file will then call the <name>_job.sh file, which will in turn call the appropriate script
  
  
## propagationMediumModels
This directory contains all files produced by P_ONE_dvirhilu/src/experimentModelCode/makeMedium.py . The folder MatthewData contains the data required to make the time residual histogram from STRAW data. The folder documentation contains details about what the medium property files contain. STRAW, ANTARES, spice_3.2.1 are medium models from STRAW, ANTARES, and IceCube respectively.

## DOMCharacteristics
contains information to describe the DOM acceptance and angular acceptance. Folders are named based on the DOM they describe.

* AngularAcceptance.dat - contains the normalized polynomial coefficients used to describe the angular acceptance. Represented using FunctionClasses.Polynomial

* DOMEfficiency.dat - first row contains a wavelength and the second row contains the DOM acceptance at that wavelenght. Values in between interpolated linearly. Represented using FunctionClasses.FunctionFromTable

## CrossSectionModels

Copied from IceCube build for convenience. Used by Neutrino-Generator
