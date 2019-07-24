#!/usr/bin/env python
'''
from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units, I3Frame
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib

infileListOrig = []
for i in range(300,399):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/HorizGeo/NuGen_step3_HorizGeo_' + str(i) + '.i3.gz')
    infileListOrig.append(infile)

infileListNew = []
for i in range(400, 499):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/HorizGeo/NuGen_step3_HorizGeo_' + str(i) + '.i3.gz')
    infileListNew.append(infile)

logEnergyOrig = []
weightsOrig = []
cosZenithOrig = []
minLogE = 0
maxLogE = 0
maxZenith = 0
minAzimuth = 0
horizontalLogEOrig = []
horizontalWeightsOrig = []
for infile in infileListOrig:
    while(infile.more()):
        frame = infile.pop_daq()
        primary = frame["NuGPrimary"]
        logEnergyOrig.append(np.log10(primary.energy))
        cosZenithOrig.append(np.cos(primary.dir.zenith))
        weightDict = frame["I3MCWeightDict"]
        oneWeight = weightDict["OneWeight"]
        # assuming all files have equal event numbers
        numEvents = weightDict["NEvents"] * len(infileListOrig)
        weightsOrig.append(oneWeight/numEvents)
        minLogE = weightDict["MinEnergyLog"]
        maxLogE = weightDict["MaxEnergyLog"]
        minAzimuth = weightDict["MinAzimuth"]
        minZenith = weightDict["MinZenith"]
        if np.cos(primary.dir.zenith) < 0.1 and np.cos(primary.dir.zenith) > -0.1:
            horizontalLogEOrig.append(np.log10(primary.energy))
            horizontalWeightsOrig.append(oneWeight/numEvents)

logEnergyNew = []
weightsNew = []
cosZenithNew = []
horizontalLogENew = []
horizontalWeightsNew = []
for infile in infileListNew:
    while(infile.more()):
        frame = infile.pop_daq()
        primary = frame["NuGPrimary"]
        logEnergyNew.append(np.log10(primary.energy))
        cosZenithNew.append(np.cos(primary.dir.zenith))
        weightDict = frame["I3MCWeightDict"]
        oneWeight = weightDict["OneWeight"]
        # assuming all files have equal event numbers
        numEvents = weightDict["NEvents"] * len(infileListNew)
        weightsNew.append(oneWeight/numEvents)
        minLogE = weightDict["MinEnergyLog"]
        maxLogE = weightDict["MaxEnergyLog"]
        minAzimuth = weightDict["MinAzimuth"]
        minZenith = weightDict["MinZenith"]
        if np.cos(primary.dir.zenith) < 0.1 and np.cos(primary.dir.zenith) > -0.1:
            horizontalLogENew.append(np.log10(primary.energy))
            horizontalWeightsNew.append(oneWeight/numEvents)

binsE = np.linspace(minLogE, maxLogE, 10)

binsZenith = np.linspace(-np.cos(minZenith),np.cos(minZenith),10)

netdE = 10**binsE[len(binsE)-1] - 10**binsE[0]

netdOmega = (binsZenith[len(binsZenith)-1] - binsZenith[0])*np.pi

dOmegaHorizon = 0.2*np.pi

dOmega = (binsZenith[1] - binsZenith[0])*np.pi
dEOrig = []
for logE in logEnergyOrig:
    position = 0
    for i in range(len(binsE)):
        if logE < binsE[i]:
            position = i
            break
    dEOrig.append(10**binsE[position] - 10**binsE[position-1])

dEHorOrig = []
for logE in horizontalLogEOrig:
    position = 0
    for i in range(len(binsE)):
        if logE < binsE[i]:
            position = i
            break
    dEHorOrig.append(10**binsE[position] - 10**binsE[position-1])

dENew = []
for logE in logEnergyNew:
    position = 0
    for i in range(len(binsE)):
        if logE < binsE[i]:
            position = i
            break
    dENew.append(10**binsE[position] - 10**binsE[position-1])

dEHorNew = []
for logE in horizontalLogENew:
    position = 0
    for i in range(len(binsE)):
        if logE < binsE[i]:
            position = i
            break
    dEHorNew.append(10**binsE[position] - 10**binsE[position-1])

areaWeights2VarOrig = [weightsOrig[i]*10**(-4)/(dEOrig[i]*dOmega) for i in range(len(weightsOrig))]
areaWeightsEnergyOrig = [weightsOrig[i]*10**(-4)/(dEOrig[i]*netdOmega) for i in range(len(weightsOrig))]
areaWeightsAngleOrig = [weightsOrig[i]*10**(-4)/(netdE*dOmega) for i in range(len(weightsOrig))]
areaWeightsHorOrig = [horizontalWeightsOrig[i]*10**(-4)/(dEHorOrig[i]*dOmegaHorizon) for i in range(len(horizontalWeightsOrig))]

areaWeights2VarNew = [weightsNew[i]*10**(-4)/(dENew[i]*dOmega) for i in range(len(weightsNew))]
areaWeightsEnergyNew = [weightsNew[i]*10**(-4)/(dENew[i]*netdOmega) for i in range(len(weightsNew))]
areaWeightsAngleNew = [weightsNew[i]*10**(-4)/(netdE*dOmega) for i in range(len(weightsNew))]
areaWeightsHorNew = [horizontalWeightsNew[i]*10**(-4)/(dEHorNew[i]*dOmegaHorizon) for i in range(len(horizontalWeightsNew))]

plt.hist(logEnergyOrig, histtype = "step", log = True, weights = areaWeightsEnergyOrig, bins = binsE)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
plt.hist(cosZenithOrig, histtype = "step", log = True, weights = areaWeightsAngleOrig, bins = binsZenith)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$\cos{\theta}$')

plt.figure()
plotOutputs = plt.hist2d(logEnergyOrig, cosZenithOrig, norm = matplotlib.colors.LogNorm(), weights = areaWeights2VarOrig, bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.colorbar(plotOutputs[3])

plt.figure()
plt.hist(horizontalLogEOrig, histtype = "step", log = True, weights = areaWeightsHorOrig, bins = binsE)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ', ' + r'$\cos{\theta}$' + ' between (-0.1,0.1) )')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.hist(logEnergyNew, histtype = "step", log = True, weights = areaWeightsEnergyNew, bins = binsE)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
plt.hist(cosZenithNew, histtype = "step", log = True, weights = areaWeightsAngleNew, bins = binsZenith)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$\cos{\theta}$')

plt.figure()
plotOutputs = plt.hist2d(logEnergyNew, cosZenithNew, norm = matplotlib.colors.LogNorm(), weights = areaWeights2VarNew, bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.colorbar(plotOutputs[3])

plt.figure()
plt.hist(horizontalLogENew, histtype = "step", log = True, weights = areaWeightsHorNew, bins = binsE)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ', ' + r'$\cos{\theta}$' + ' between (-0.1,0.1) )')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
h1, xedges, yedges = np.histogram2d(logEnergyOrig, cosZenithOrig, bins = [binsE, binsZenith], weights = areaWeights2VarOrig)
h2, xedges, yedges = np.histogram2d(logEnergyNew, cosZenithNew, bins = [binsE, binsZenith], weights = areaWeights2VarNew)
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
plt.title("Effective Area Ratio")
plt.show()
'''

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib

infileListOrig = []
for i in range(300,399):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/HorizGeo/NuGen_step3_HorizGeo_' + str(i) + '.i3.gz')
    infileListOrig.append(infile)

logEnergyOrig = []
weightsOrig = []
cosZenithOrig = []
minLogE = 0
maxLogE = 0
maxZenith = 0
minAzimuth = 0
horizontalLogEOrig = []
horizontalWeightsOrig = []
logEnergyNew = []
weightsNew = []
cosZenithNew = []
horizontalLogENew = []
horizontalWeightsNew = []
for infile in infileListOrig:
    while(infile.more()):
        frame = infile.pop_daq()
        primary = frame["NuGPrimary"]
        logEnergyOrig.append(np.log10(primary.energy))
        cosZenithOrig.append(np.cos(primary.dir.zenith))
        weightDict = frame["I3MCWeightDict"]
        oneWeight = weightDict["OneWeight"]
        # assuming all files have equal event numbers
        numEvents = weightDict["NEvents"] * len(infileListOrig)
        weightsOrig.append(oneWeight/numEvents)
        minLogE = weightDict["MinEnergyLog"]
        maxLogE = weightDict["MaxEnergyLog"]
        minAzimuth = weightDict["MinAzimuth"]
        minZenith = weightDict["MinZenith"]
        if np.cos(primary.dir.zenith) < 0.1 and np.cos(primary.dir.zenith) > -0.1:
            horizontalLogEOrig.append(np.log10(primary.energy))
            horizontalWeightsOrig.append(oneWeight/numEvents)
        mcpeMap = frame["MCPESeriesMap"]
        domCount = 0
        for mcpeList in mcpeMap.values():
            if len(mcpeList) >= 6:
                domCount += 1
    
        if domCount >= 6:
            logEnergyNew.append(np.log10(primary.energy))
            cosZenithNew.append(np.cos(primary.dir.zenith))
            weightsNew.append(oneWeight/numEvents)
            if np.cos(primary.dir.zenith) < 0.1 and np.cos(primary.dir.zenith) > -0.1:
                horizontalLogENew.append(np.log10(primary.energy))
                horizontalWeightsNew.append(oneWeight/numEvents)


binsE = np.linspace(minLogE, maxLogE, 10)

binsZenith = np.linspace(-np.cos(minZenith),np.cos(minZenith),10)

netdE = 10**binsE[len(binsE)-1] - 10**binsE[0]

netdOmega = (binsZenith[len(binsZenith)-1] - binsZenith[0])*np.pi

dOmegaHorizon = 0.2*np.pi

dOmega = (binsZenith[1] - binsZenith[0])*np.pi
dEOrig = []
for logE in logEnergyOrig:
    position = 0
    for i in range(len(binsE)):
        if logE < binsE[i]:
            position = i
            break
    dEOrig.append(10**binsE[position] - 10**binsE[position-1])

dEHorOrig = []
for logE in horizontalLogEOrig:
    position = 0
    for i in range(len(binsE)):
        if logE < binsE[i]:
            position = i
            break
    dEHorOrig.append(10**binsE[position] - 10**binsE[position-1])

dENew = []
for logE in logEnergyNew:
    position = 0
    for i in range(len(binsE)):
        if logE < binsE[i]:
            position = i
            break
    dENew.append(10**binsE[position] - 10**binsE[position-1])

dEHorNew = []
for logE in horizontalLogENew:
    position = 0
    for i in range(len(binsE)):
        if logE < binsE[i]:
            position = i
            break
    dEHorNew.append(10**binsE[position] - 10**binsE[position-1])

print(len(logEnergyNew), len(logEnergyOrig))

areaWeights2VarOrig = [weightsOrig[i]*10**(-4)/(dEOrig[i]*dOmega) for i in range(len(weightsOrig))]
areaWeightsEnergyOrig = [weightsOrig[i]*10**(-4)/(dEOrig[i]*netdOmega) for i in range(len(weightsOrig))]
areaWeightsAngleOrig = [weightsOrig[i]*10**(-4)/(netdE*dOmega) for i in range(len(weightsOrig))]
areaWeightsHorOrig = [horizontalWeightsOrig[i]*10**(-4)/(dEHorOrig[i]*dOmegaHorizon) for i in range(len(horizontalWeightsOrig))]

areaWeights2VarNew = [weightsNew[i]*10**(-4)/(dENew[i]*dOmega) for i in range(len(weightsNew))]
areaWeightsEnergyNew = [weightsNew[i]*10**(-4)/(dENew[i]*netdOmega) for i in range(len(weightsNew))]
areaWeightsAngleNew = [weightsNew[i]*10**(-4)/(netdE*dOmega) for i in range(len(weightsNew))]
areaWeightsHorNew = [horizontalWeightsNew[i]*10**(-4)/(dEHorNew[i]*dOmegaHorizon) for i in range(len(horizontalWeightsNew))]
'''
plt.hist(logEnergyOrig, histtype = "step", log = True, weights = areaWeightsEnergyOrig, bins = binsE)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
plt.hist(cosZenithOrig, histtype = "step", log = True, weights = areaWeightsAngleOrig, bins = binsZenith)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$\cos{\theta}$')

plt.figure()
plotOutputs = plt.hist2d(logEnergyOrig, cosZenithOrig, norm = matplotlib.colors.LogNorm(), weights = areaWeights2VarOrig, bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.colorbar(plotOutputs[3])

plt.figure()
plt.hist(horizontalLogEOrig, histtype = "step", log = True, weights = areaWeightsHorOrig, bins = binsE)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ', ' + r'$\cos{\theta}$' + ' between (-0.1,0.1) )')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
plt.hist(logEnergyNew, histtype = "step", log = True, weights = areaWeightsEnergyNew, bins = binsE)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
plt.hist(cosZenithNew, histtype = "step", log = True, weights = areaWeightsAngleNew, bins = binsZenith)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.xlabel(r'$\cos{\theta}$')

plt.figure()
plotOutputs = plt.hist2d(logEnergyNew, cosZenithNew, norm = matplotlib.colors.LogNorm(), weights = areaWeights2VarNew, bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ')')
plt.colorbar(plotOutputs[3])

plt.figure()
plt.hist(horizontalLogENew, histtype = "step", log = True, weights = areaWeightsHorNew, bins = binsE)
plt.title("Effective Area Distribution (in " + r'$m^2$' + ', ' + r'$\cos{\theta}$' + ' between (-0.1,0.1) )')
plt.xlabel(r'$log_{10}\, E/GeV$')

plt.figure()
h1, xedges, yedges = np.histogram2d(logEnergyOrig, cosZenithOrig, bins = [binsE, binsZenith], weights = areaWeights2VarOrig)
h2, xedges, yedges = np.histogram2d(logEnergyNew, cosZenithNew, bins = [binsE, binsZenith], weights = areaWeights2VarNew)
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
plt.title("Effective Area Ratio")
'''
plt.figure()
h1, edges= np.histogram(logEnergyNew, bins = binsE, weights = areaWeightsEnergyNew)
h2, edges = np.histogram(logEnergyOrig, bins = binsE, weights = areaWeightsEnergyOrig)
for i in range(len(edges)-1):
    if h2[i] <=0.0000001:
        h2[i] = 1.0
        h1[i] = 0
print h1
print h2
print edges
ratioHist = h1 / h2
plt.bar(edges[:-1], ratioHist, align = 'edge')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$ratio$')
plt.title("Ratio Of Effective Areas (5 DOMs, 3 Hits/DOM)")

plt.figure()
h1, edges= np.histogram(horizontalLogENew, bins = binsE, weights = areaWeightsHorNew)
h2, edges = np.histogram(horizontalLogEOrig, bins = binsE, weights = areaWeightsHorOrig)
for i in range(len(edges)-1):
    if h2[i] <=0.0000001:
        h2[i] = 1.0
        h1[i] = 0
print h1
print h2
print edges
print h1/h2
ratioHist = h1 / h2
plt.bar(edges[:-1], ratioHist, align = 'edge')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$ratio$')
plt.title("Ratio Of Effective Areas, cos zenith between (-0.1,0.1)")

plt.show()
