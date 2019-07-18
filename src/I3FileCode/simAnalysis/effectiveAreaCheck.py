#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units, I3Frame
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib

infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/examples/Level2_Pass3_IC86.2016_NuMu.020808.000280.i3.zst')

logEnergy = []
cosZenith = []
weights = []
minLogE = 0
maxLogE = 0

for frame in infile:
    if frame.Stop == I3Frame.DAQ:
        primary = frame["I3MCTree_preMuonProp"].primaries[0]
        logEnergy.append(np.log10(primary.energy))
        cosZenith.append(np.cos(primary.dir.zenith))
        weightDict = frame["I3MCWeightDict"]
        oneWeight = weightDict["OneWeight"]
        # assuming all files have equal event numbers
        numEvents = weightDict["NEvents"]
        weights.append(oneWeight/numEvents)
        minLogE = weightDict["MinEnergyLog"]
        maxLogE = weightDict["MaxEnergyLog"]


binsE = np.linspace(minLogE, maxLogE, 10)
binsZenith = np.linspace(-1, 1, 10)

dOmega = (binsZenith[1] - binsZenith[0])*np.pi*2
dE = []
for logE in logEnergy:
    position = 0
    for i in range(len(binsE)):
        if logE < binsE[i]:
            position = i
            break
    dE.append(10**binsE[position] - 10**binsE[position-1])

areaWeights = [weights[i]*10**(-9)/(dE[i]*dOmega) for i in range(len(weights))]

plt.hist(logEnergy, histtype = "step", log = True, weights = areaWeights, bins = binsE)
plt.title("Effective Area Distribution (in " + r'$cm^2$' + ')')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
plt.hist(cosZenith, histtype = "step", log = True, weights = areaWeights, bins = binsZenith)
plt.title("Effective Area Distribution (in " + r'$cm^2$' + ')')
plt.xlabel(r'$\cos{\theta}$')

plt.figure()
plotOutputs = plt.hist2d(logEnergy,cosZenith, norm = matplotlib.colors.LogNorm(), weights = areaWeights, bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$cm^2$' + ')')
plt.colorbar(plotOutputs[3])

plt.show()