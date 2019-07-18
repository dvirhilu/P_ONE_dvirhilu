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

infileListStep2 = []
for i in range(int(args.minFileNum), int(args.maxFileNum) + 1):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep2/NuGen_step2_' + str(args.runType) + '_' + str(i) + '.i3.gz')
    infileListStep2.append(infile)


infileListStep3 = []
for i in range(int(args.minFileNum), int(args.maxFileNum) + 1):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/NuGen_step3_' + str(args.runType) + '_' + str(i) + '.i3.gz')
    infileListStep3.append(infile)

energyStep1 = []
zenithStep1 = []
weights = []
for infile in infileListStep1:
    while infile.more():
        frame = infile.pop_daq()
        primary = frame["NuGPrimary"]
        weights.append(frame["EventWeight"].value/len(infileListStep1))
        zenithStep1.append(primary.dir.zenith)
        energyStep1.append(primary.energy)

logEStep1 = np.log10(energyStep1)
cosZenStep1 = np.cos(zenithStep1)

'''
plt.hist(logEStep1, histtype = "step", log = True, weights = weights, bins = 30)
plt.title("Weighted Neutrino Energy Distribution")
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
plt.hist(cosZenStep1, histtype = "step", log = True, weights = weights, bins = 30)
plt.title("Weighted Muon Angular Distribution (Zenith)")
plt.xlabel(r'$\cos{\theta}$')

plt.figure()
plotOutputs = plt.hist2d(logEStep1,cosZenStep1, norm = matplotlib.colors.LogNorm(), bins = 10)
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Neutrino Energy and Zenith Distribution")
plt.colorbar(plotOutputs[3])

plt.show()
'''

logEStep2 = []
cosZenStep2 = []
for infile in infileListStep2:
    while infile.more():
        frame = infile.pop_daq()
        primary = frame["NuGPrimary"]
        cosZenStep2.append(np.cos(primary.dir.zenith))
        logEStep2.append(np.log10(primary.energy))


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

areaWeights = [weights[i]*10**(-4)/(dE[i]*dOmega) for i in range(len(weights))]
print binsE
print binsZenith
print areaWeights

plt.hist(logEnergy, histtype = "step", log = True, weights = areaWeights, bins = binsE)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
plt.hist(cosZenith, histtype = "step", log = True, weights = areaWeights, bins = binsZenith)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$\cos{\theta}$')

plt.figure()
plotOutputs = plt.hist2d(logEnergy,cosZenith, norm = matplotlib.colors.LogNorm(), weights = areaWeights, bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.colorbar(plotOutputs[3])

#plt.show()

plt.figure()
plotouts = plt.hist2d(logEStep1, cosZenStep1, norm = matplotlib.colors.LogNorm(), bins = [binsE, binsZenith])
plt.colorbar(plotouts[3])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Event Distribution - Step 1")

plt.figure()
plotouts = plt.hist2d(logEStep2, cosZenStep2, norm = matplotlib.colors.LogNorm(), bins = [binsE, binsZenith])
plt.colorbar(plotouts[3])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Event Distribution - Step 2")

plt.figure()
plotouts = plt.hist2d(logEnergy, cosZenith, norm = matplotlib.colors.LogNorm(), bins = [binsE, binsZenith])
plt.colorbar(plotouts[3])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Event Distribution - Step 3")

plt.figure()
h1, xedges, yedges = np.histogram2d(logEStep1, cosZenStep1, bins = [binsE, binsZenith])
h2, xedges, yedges = np.histogram2d(logEnergy, cosZenith, bins = [binsE, binsZenith])
for i in range(len(xedges)-1):
    for j in range(len(yedges)-1):
        if h1[j][i] <=0.0000001:
            h1[j][i] = 1.0
            h2[j][i] = 0
ratioHist = h2 / h1
pc = plt.pcolor(xedges, yedges, ratioHist.T, norm = matplotlib.colors.LogNorm())
plt.colorbar(pc)
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Detection Efficiency")
plt.show()
