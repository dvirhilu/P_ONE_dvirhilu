#!/usr/bin/env python

from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame
from icecube.dataclasses import I3Particle
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "Creates a reconstruction of the muon track using a linear least squares fit on the pulses")
parser.add_argument( '-i', '--infile', help = "input file used" )
args = parser.parse_args()

infile = dataio.I3File(str(args.infile))

cosAlpha = []
logEnergy = []
#i = 0
for frame in infile:
    primary = frame["NuGPrimary"]
    mctree = frame["I3MCTree"]
    muon = dataclasses.I3MCTree.first_child(mctree, primary)
    recoParticle = frame["LineFitRecoParticle"]

    muonDir = muon.dir
    recoDir = recoParticle.dir

    dotProduct = muonDir.x*recoDir.x + muonDir.y*recoDir.y + muonDir.z*recoDir.z
    cosAlpha.append(dotProduct)

    logEnergy.append(np.log10(primary.energy))    

    #if dotProduct < 0.1:
    #    print i
    #i += 1


alpha = [np.arccos(cosA)/I3Units.deg for cosA in cosAlpha]

plt.hist(cosAlpha, histtype = 'step', bins = 20)
plt.xlabel(r'$\cos{\alpha}$')
plt.title('Distribution of Relative Angle of Muon and its Reconstruction')

plt.figure()
plt.hist(cosAlpha, histtype = 'step', bins = 20, log = True)
plt.xlabel(r'$\cos{\alpha}$')
plt.title('Distribution of Relative Angle of Muon and its Reconstruction')

plt.figure()
plt.hist(alpha, histtype = 'step', bins = 20)
plt.xlabel(r'$\alpha$')
plt.title('Distribution of Relative Angle of Muon and its Reconstruction')

plt.figure()
plt.hist(alpha, histtype = 'step', bins = 20, log = True)
plt.xlabel(r'$\alpha}$')
plt.title('Distribution of Relative Angle of Muon and its Reconstruction')

plt.figure()
h, xedges, yedges = np.histogram2d(logEnergy, alpha)
for i in range(len(xedges)-1):
    total = 0
    for j in range(len(yedges)-1):
        total += h[i][j]
    if total == 0:
        h[i] = [0*element for element in h[i]]
    else:
        h[i] = [element/total for element in h[i]]

pc = plt.pcolor(xedges, yedges, h.T)
plt.colorbar(pc)
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel('Relative Angle')
plt.title("Detection Efficiency")
plt.show()