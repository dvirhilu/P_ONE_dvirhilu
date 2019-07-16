#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib

parser = argparse.ArgumentParser(description = "Find neutrino distribution from NuGen simulation")
parser.add_argument( '-n', '--minFileNum', help = "smallest file number used" )
parser.add_argument( '-N', '--maxFileNum', help = "largest file number used")
parser.add_argument( '-r', '--runType', help = "The run type (must align with GCD type)")
args = parser.parse_args()

# open file
infileListStep1 = []
for i in range(int(args.minFileNum), int(args.maxFileNum) + 1):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep1/NuGen_step1_' + str(args.runType) + '_' + str(i) + '.i3.gz')
    infileListStep1.append(infile)

infileListStep3 = []
for i in range(int(args.minFileNum), int(args.maxFileNum) + 1):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/NuGen_step3_' + str(args.runType) + '_' + str(i) + '.i3.gz')
    infileListStep3.append(infile)

energyStep1 = []
zenithStep1 = []
azimuthStep1 = []
weights = []
for infile in infileListStep1:
    while infile.more():
        frame = infile.pop_daq()
        primary = frame["NuGPrimary"]
        weights.append(frame["EventWeight"].value/len(infileListStep1))
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


logEnergy = []
cosZenith = []
weights = []
minLogE = 0
maxLogE = 0
for infile in infileListStep3:
    while(infile.more()):
        frame = infile.pop_daq()
        primary = frame["NuGPrimary"]
        logEnergy.append(np.log10(primary.energy))
        cosZenith.append(np.cos(primary.dir.zenith))
        weightDict = frame["I3MCWeightDict"]
        oneWeight = weightDict["OneWeight"]
        # assuming all files have equal event numbers
        numEvents = weightDict["NEvents"] * len(infileListStep3)
        weights.append(oneWeight/numEvents/2)
        minLogE = weightDict["MinEnergyLog"]
        maxLogE = weightDict["MaxEnergyLog"]

binsE = np.logspace(minLogE, maxLogE, 10, base = 10.0)
binsZenith = np.linspace(-1, 1, 10)
print(numEvents)

dOmega = (binsZenith[1] - binsZenith[0])*np.pi*2
dE = []
for logE in logEnergy:
    position = 0
    for i in range(len(binsE)):
        if logE < np.log10(binsE[i]):
            print(logE)
            print(np.log10(binsE[i]))
            position = i
            break
    dE.append(binsE[position] - binsE[position-1])

areaWeights = [weights[i]*10**(-9)/(dE[i]*dOmega) for i in range(len(weights))]

plt.hist(logEnergy, histtype = "step", log = True, weights = areaWeights)
plt.title("Effective Area Distribution (in " + r'$cm^2$' + ')')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
plt.hist(cosZenith, histtype = "step", log = True, weights = areaWeights)
plt.title("Effective Area Distribution (in " + r'$cm^2$' + ')')
plt.xlabel(r'$\cos{\theta}$')

plt.figure()
plotOutputs = plt.hist2d(logEnergy,cosZenith, norm = matplotlib.colors.LogNorm(), weights = areaWeights)
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$cm^2$' + ')')
plt.colorbar(plotOutputs[3])

plt.show()
