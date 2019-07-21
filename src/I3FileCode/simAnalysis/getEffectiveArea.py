#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units, I3Frame
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
cosRelAnglStep1 = []
for infile in infileListStep1:
    while infile.more():
        frame = infile.pop_daq()
        primary = frame["NuGPrimary"]
        cosRelAnglStep1.append(-primary.dir.x)
        energyStep1.append(primary.energy)

logEStep1 = np.log10(energyStep1)

logEStep2 = []
cosRelAnglStep2 = []
for infile in infileListStep2:
    while infile.more():
        frame = infile.pop_daq()
        primary = frame["NuGPrimary"]
        cosRelAnglStep2.append(-primary.dir.x)
        logEStep2.append(np.log10(primary.energy))


logEnergy = []
cosRelAngl = []
weights = []
minLogE = 0
maxLogE = 0
maxZenith = 0
minAzimuth = 0
for infile in infileListStep3:
    while(infile.more()):
        frame = infile.pop_daq()
        primary = frame["NuGPrimary"]
        logEnergy.append(np.log10(primary.energy))
        cosRelAngl.append(-primary.dir.x)
        weightDict = frame["I3MCWeightDict"]
        oneWeight = weightDict["OneWeight"]
        # assuming all files have equal event numbers
        numEvents = weightDict["NEvents"] * len(infileListStep3)
        weights.append(oneWeight/numEvents)
        minLogE = weightDict["MinEnergyLog"]
        maxLogE = weightDict["MaxEnergyLog"]
        minAzimuth = weightDict["MinAzimuth"]
        minZenith = weightDict["MinZenith"]

print len(logEnergy)

binsE = np.linspace(minLogE, maxLogE, 10)
binsRelAngl = np.linspace(np.sin(minZenith)*np.cos(minAzimuth), 1, 10)

netdE = 10**binsE[len(binsE)-1] - 10**binsE[0]
netdOmega = binsRelAngl[len(binsRelAngl)-1] - binsRelAngl[0]

dOmega = (binsRelAngl[1] - binsRelAngl[0])*np.pi*2
dE = []
for logE in logEnergy:
    position = 0
    for i in range(len(binsE)):
        if logE < binsE[i]:
            position = i
            break
    dE.append(10**binsE[position] - 10**binsE[position-1])

areaWeights2Var = [weights[i]*10**(-4)/(dE[i]*dOmega) for i in range(len(weights))]
areaWeightsEnergy = [weights[i]*10**(-4)/(dE[i]*netdOmega) for i in range(len(weights))]
areaWeightsAngle = [weights[i]*10**(-4)/(netdE*dOmega) for i in range(len(weights))]

plt.hist(logEnergy, histtype = "step", log = True, weights = areaWeightsEnergy, bins = binsE)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
plt.hist(cosRelAngl, histtype = "step", log = True, weights = areaWeightsAngle, bins = binsRelAngl)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$\cos{\theta_{rel}}$')

plt.figure()
plotOutputs = plt.hist2d(logEnergy,cosRelAngl, norm = matplotlib.colors.LogNorm(), weights = areaWeights2Var, bins = [binsE, binsRelAngl])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta_{rel}}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.colorbar(plotOutputs[3])

#plt.show()

plt.figure()
plotouts = plt.hist2d(logEStep1, cosRelAnglStep1, norm = matplotlib.colors.LogNorm(), bins = [binsE, binsRelAngl])
plt.colorbar(plotouts[3])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta_{rel}}$')
plt.title("Event Distribution - Step 1")

plt.figure()
plotouts = plt.hist2d(logEStep2, cosRelAnglStep2, norm = matplotlib.colors.LogNorm(), bins = [binsE, binsRelAngl])
plt.colorbar(plotouts[3])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta_{rel}}$')
plt.title("Event Distribution - Step 2")

plt.figure()
plotouts = plt.hist2d(logEnergy, cosRelAngl, norm = matplotlib.colors.LogNorm(), bins = [binsE, binsRelAngl])
plt.colorbar(plotouts[3])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta_{rel}}$')
plt.title("Event Distribution - Step 3")

plt.figure()
h1, xedges, yedges = np.histogram2d(logEStep1, cosRelAnglStep1, bins = [binsE, binsRelAngl])
h2, xedges, yedges = np.histogram2d(logEnergy, cosRelAngl, bins = [binsE, binsRelAngl])
for i in range(len(xedges)-1):
    for j in range(len(yedges)-1):
        if h1[j][i] <=0.0000001:
            h1[j][i] = 1.0
            h2[j][i] = 0
ratioHist = h2 / h1
pc = plt.pcolor(xedges, yedges, ratioHist.T, norm = matplotlib.colors.LogNorm())
plt.colorbar(pc)
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta_{rel}}$')
plt.title("Detection Efficiency")
plt.show()

infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/examples/Level2_Pass3_IC86.2016_NuMu.020808.000280.i3.zst')

logEnergyIceCube = []
cosRelAnglIceCube = []
weightsIceCube = []
minLogE = 0
maxLogE = 0

for frame in infile:
    if frame.Stop == I3Frame.DAQ:
        primary = frame["I3MCTree_preMuonProp"].primaries[0]
        logEnergyIceCube.append(np.log10(primary.energy))
        cosRelAnglIceCube.append(primary.dir.x)
        weightDict = frame["I3MCWeightDict"]
        oneWeight = weightDict["OneWeight"]
        # assuming all files have equal event numbers
        numEvents = weightDict["NEvents"]
        weightsIceCube.append(oneWeight/numEvents)
        minLogE = weightDict["MinEnergyLog"]
        maxLogE = weightDict["MaxEnergyLog"]

print len(logEnergyIceCube)

dE = []
for logE in logEnergyIceCube:
    position = 0
    for i in range(len(binsE)):
        if logE < binsE[i]:
            position = i
            break
    dE.append(10**binsE[position] - 10**binsE[position-1])

areaWeights2VarIC = [weightsIceCube[i]*10**(-4)/(dE[i]*dOmega) for i in range(len(weightsIceCube))]
areaWeightsEnergyIC = [weightsIceCube[i]*10**(-4)/(dE[i]*netdOmega) for i in range(len(weightsIceCube))]
areaWeightsAngleIC = [weightsIceCube[i]*10**(-4)/(netdE*dOmega) for i in range(len(weightsIceCube))]

plt.hist(logEnergyIceCube, histtype = "step", log = True, weights = areaWeightsEnergyIC, bins = binsE)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
plt.hist(cosRelAnglIceCube, histtype = "step", log = True, weights = areaWeightsAngleIC, bins = binsRelAngl)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$\cos{\theta_{rel}}$')

plt.figure()
plotOutputs = plt.hist2d(logEnergyIceCube,cosRelAnglIceCube, norm = matplotlib.colors.LogNorm(), weights = areaWeights2VarIC, bins = [binsE, binsRelAngl])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta_{rel}}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.colorbar(plotOutputs[3])

plt.figure()
h1, edges= np.histogram(logEnergy, bins = binsE, weights = areaWeightsEnergy)
h2, edges = np.histogram(logEnergyIceCube, bins = binsE, weights = areaWeightsEnergyIC)
for i in range(len(edges)-1):
    if h2[i] <=0.0000001:
        h2[i] = 1.0
        h1[i] = 0
print h1
print h2
print edges
ratioHist = h1 / h2
plt.bar(edges[:-1], ratioHist, log = True, align = 'edge')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$ratio$')
plt.title("Ratio Of Effective Areas")
plt.show()