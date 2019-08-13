#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib
from simAnalysis import SimAnalysis

parser = argparse.ArgumentParser(description = "Creates a reconstruction of the muon track using a linear least squares fit on the pulses")
parser.add_argument( '-H', '--hitThresh', default = 6, help = "threshold of hits for the DOM to be considered")
parser.add_argument( '-D', '--domThresh', default = 6 , help = "threshold of hit DOMs for the frame to be considered")
args = parser.parse_args()

hitThresh = int(args.hitThresh)
domThresh = int(args.domThresh)

infileListOrig = []
for i in range(300, 400):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/HorizGeo/NuGen_step3_HorizGeo_' + str(i) + '.i3.gz')
    infileListOrig.append(infile)

gcdOrig = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/corHorizgeo/CorrHorizGeo_n15_b100.0_a18.0_l3_rise_fall_offset_simple_spacing.i3.gz')
domsUsedOrig = gcdOrig.pop_frame()["I3Geometry"].omgeo.keys()

infileListNew = []
for i in range(500, 700):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/denseGeo/NuGen_step3_denseGeo_' + str(i) + '.i3.gz')
    infileListNew.append(infile)

gcdl0 = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/partialDenseGeo/partialDenseGeo_comparisonGeometry_l0.i3.gz')
gcdl1 = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/partialDenseGeo/partialDenseGeo_comparisonGeometry_l1.i3.gz')
gcdl2 = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/partialDenseGeo/partialDenseGeo_comparisonGeometry_l2.i3.gz')
domsUsedl0 = gcdl0.pop_frame()["I3Geometry"].omgeo.keys()
domsUsedl1 = gcdl1.pop_frame()["I3Geometry"].omgeo.keys()
domsUsedl2 = gcdl2.pop_frame()["I3Geometry"].omgeo.keys()

dataOrig, weightsOrig, binsOrig = SimAnalysis.getEffectiveAreaData(infileListOrig, domsUsedOrig, hitThresh, domThresh, 10)
print "finished orig"
datal0, weightsl0, _binsl0 = SimAnalysis.getEffectiveAreaData(infileListNew, domsUsedl0, hitThresh, domThresh, 10)
print "finished l0"
for infile in infileListNew:
    infile.rewind()
datal1, weightsl1, _binsl1 = SimAnalysis.getEffectiveAreaData(infileListNew, domsUsedl1, hitThresh, domThresh, 10)
print "finished l1"
for infile in infileListNew:
    infile.rewind()
datal2, weightsl2, _binsl2 = SimAnalysis.getEffectiveAreaData(infileListNew, domsUsedl2, hitThresh, domThresh, 10)
print "finished l2"

binsE = binsOrig["logEnergy"]
binsZenith = binsOrig["cosZenith"]

plt.figure()
plotOutputs = plt.hist2d(dataOrig["logEnergy"], dataOrig["cosZenith"], norm = matplotlib.colors.LogNorm(), weights = weightsOrig["2Variable"], bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ') Original')
plt.colorbar(plotOutputs[3])

plt.figure()
plotOutputs = plt.hist2d(datal0["logEnergy"], datal0["cosZenith"], norm = matplotlib.colors.LogNorm(), weights = weightsl0["2Variable"], bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ') l0')
plt.colorbar(plotOutputs[3])

plt.figure()
plotOutputs = plt.hist2d(datal1["logEnergy"], datal1["cosZenith"], norm = matplotlib.colors.LogNorm(), weights = weightsl1["2Variable"], bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ') l1')
plt.colorbar(plotOutputs[3])

plt.figure()
plotOutputs = plt.hist2d(datal2["logEnergy"], datal2["cosZenith"], norm = matplotlib.colors.LogNorm(), weights = weightsl2["2Variable"], bins = [binsE, binsZenith])
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$\cos{\theta}$')
plt.title("Effective Area Distribution (in " + r'$m^2$' + ') l2')
plt.colorbar(plotOutputs[3])


'''
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
h1, edges= np.histogram(dataOrig["logEnergy"], bins = binsE, weights = weightsOrig["logEnergy"])
h2, edges = np.histogram(datal0["logEnergy"], bins = binsE, weights = weightsl0["logEnergy"])
for i in range(len(edges)-1):
    if h2[i] <=0.0000001:
        h2[i] = 1.0
        h1[i] = 0

ratioHist = h1 / h2
plt.step(edges[:-1], ratioHist, where = 'post')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$ratio$')
plt.title("Ratio Of orig/l0 (" + str(args.domThresh) + " DOMs, " + str(args.hitThresh) + " Hits/DOM)")

plt.figure()
h1, edges= np.histogram(dataOrig["Horizontal"], bins = binsE, weights = weightsOrig["Horizontal"])
h2, edges = np.histogram(datal0["Horizontal"], bins = binsE, weights = weightsl0["Horizontal"])
for i in range(len(edges)-1):
    if h2[i] <=0.0000001:
        h2[i] = 1.0
        h1[i] = 0

ratioHist = h1 / h2
plt.step(edges[:-1], ratioHist, where = 'post')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$ratio$')
plt.title("Ratio Of orig/l0 (" + str(args.domThresh) + " DOMs, " + str(args.hitThresh) + " Hits/DOM), " + r'$cos_{\theta} \epsilon(-0.1,0.1)$')



plt.figure()
h1, edges= np.histogram(dataOrig["logEnergy"], bins = binsE, weights = weightsOrig["logEnergy"])
h2, edges = np.histogram(datal1["logEnergy"], bins = binsE, weights = weightsl1["logEnergy"])
for i in range(len(edges)-1):
    if h2[i] <=0.0000001:
        h2[i] = 1.0
        h1[i] = 0

ratioHist = h1 / h2
plt.step(edges[:-1], ratioHist, where = 'post')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$ratio$')
plt.title("Ratio Of orig/l1 (" + str(args.domThresh) + " DOMs, " + str(args.hitThresh) + " Hits/DOM)")

plt.figure()
h1, edges= np.histogram(dataOrig["Horizontal"], bins = binsE, weights = weightsOrig["Horizontal"])
h2, edges = np.histogram(datal1["Horizontal"], bins = binsE, weights = weightsl1["Horizontal"])
for i in range(len(edges)-1):
    if h2[i] <=0.0000001:
        h2[i] = 1.0
        h1[i] = 0

ratioHist = h1 / h2
plt.step(edges[:-1], ratioHist, where = 'post')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$ratio$')
plt.title("Ratio Of orig/l1 (" + str(args.domThresh) + " DOMs, " + str(args.hitThresh) + " Hits/DOM), " + r'$cos_{\theta} \epsilon(-0.1,0.1)$')



plt.figure()
h1, edges= np.histogram(dataOrig["logEnergy"], bins = binsE, weights = weightsOrig["logEnergy"])
h2, edges = np.histogram(datal2["logEnergy"], bins = binsE, weights = weightsl2["logEnergy"])
for i in range(len(edges)-1):
    if h2[i] <=0.0000001:
        h2[i] = 1.0
        h1[i] = 0

ratioHist = h1 / h2
plt.step(edges[:-1], ratioHist, where = 'post')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$ratio$')
plt.title("Ratio Of orig/l2 (" + str(args.domThresh) + " DOMs, " + str(args.hitThresh) + " Hits/DOM)")

plt.figure()
h1, edges= np.histogram(dataOrig["Horizontal"], bins = binsE, weights = weightsOrig["Horizontal"])
h2, edges = np.histogram(datal2["Horizontal"], bins = binsE, weights = weightsl2["Horizontal"])
for i in range(len(edges)-1):
    if h2[i] <=0.0000001:
        h2[i] = 1.0
        h1[i] = 0

ratioHist = h1 / h2
plt.step(edges[:-1], ratioHist, where = 'post')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel(r'$ratio$')
plt.title("Ratio Of orig/l2 (" + str(args.domThresh) + " DOMs, " + str(args.hitThresh) + " Hits/DOM), " + r'$cos_{\theta} \epsilon(-0.1,0.1)$')

plt.show()
