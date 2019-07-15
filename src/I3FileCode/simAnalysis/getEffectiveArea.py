#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib

parser = argparse.ArgumentParser(description = "Find neutrino distribution from NuGen simulation")
parser.add_argument( '-n', '--minFileNum', help = "smallest file number used" )
parser.add_argument( '-N', '--maxFileNum', help = "largest file number used")
parser.add_argument( '-r', '--runType', help = "The run type (must align with sim run type)")
args = parser.parse_args()

# open file
infileListStep1 = []
for i in range(int(args.minFileNum), int(args.maxFileNum) + 1):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep1/NuGen_step1_' + str(args.runType) + '_' + str(i) + '.i3.gz')
    infileListStep1.append(infile)
    print('file ' + str(i) + 'done' )


energyStep1 = []
zenithStep1 = []
azimuthStep1 = []
weights = []
for infile in infileListStep1:
    while infile.more():
        frame = infile.pop_daq()
        primary = frame["NuGPrimary"]
        weights.append(frame["EventWeight"].value)
        zenithStep1.append(primary.dir.zenith)
        azimuthStep1.append(primary.dir.azimuth)
        energyStep1.append(primary.energy)


logEStep1 = np.log10(energyStep1)
cosZenStep1 = np.cos(zenithStep1)
azimuthStep1 = [angle/I3Units.deg for angle in azimuthStep1]

plt.hist(logEStep1, histtype = "step", log = True, weights = weights, bins = 30)
plt.title("Weighted Neutrino Energy Distribution")
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
plt.hist(azimuthStep1, histtype = "step", log = True, weights = weights, bins = 30)
plt.title("Weighted Muon Angular Distribution (Azimuth)")
plt.xlabel("Azimuth Angle (degrees)")

plt.figure()
plt.hist(cosZenStep1, histtype = "step", log = True, weights = weights, bins = 30)
plt.title("Weighted Muon Angular Distribution (Zenith)")
plt.xlabel("Cosine of the Zenith Angle")

plt.figure()
plotOutputs = plt.hist2d(logEStep1,cosZenStep1, norm = matplotlib.colors.LogNorm(), bins = 10)
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Neutrino Energy and Zenith Distribution")
plt.colorbar(plotOutputs[3])

plt.show()

