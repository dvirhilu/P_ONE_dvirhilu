#!/usr/bin/env python

from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame
from icecube.dataclasses import I3Particle
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


infileImproved = dataio.I3File('/home/dvir/workFolder/I3Files/improvedReco/HorizGeo/NuGen_improvedReco_HorizGeo_improvedRecoTest.i3.gz')
infileLineFit = dataio.I3File('/home/dvir/workFolder/I3Files/linefitReco/HorizGeo/NuGen_linefitReco_HorizGeo_testNewAlgorithm.i3.gz')

def findIndex(energy, bins):
    logE = np.log10(energy)

    for i in range(len(bins)):
        if logE < bins[i+1]:
            return i
    
    raise ValueError("energy not in range")


def fitanalysis(infile, binsE, fitType):
    cosAlpha = []
    listsError = [[] for i in range(10)]
    unsuccessfulRecoEnergy = []
    successfulRecoEnergy = []

    for frame in infile:
        primary = frame["NuGPrimary"]
        mctree = frame["I3MCTree"]
        muon = dataclasses.I3MCTree.first_child(mctree, primary)
        
        if fitType == "linefit":
            recoParticle = frame["LineFitRecoParticle"]
        else:
            recoParticle = frame["ImprovedRecoParticle"]

        if recoParticle.fit_status == dataclasses.I3Particle.InsufficientQuality:
            unsuccessfulRecoEnergy.append(np.log10(primary.energy))
            continue

        successfulRecoEnergy.append(np.log10(primary.energy))
        muonDir = muon.dir
        recoDir = recoParticle.dir

        dotProduct = muonDir.x*recoDir.x + muonDir.y*recoDir.y + muonDir.z*recoDir.z
        cosAlpha.append(dotProduct)

        index = findIndex(primary.energy, binsE)
        listsError[index].append(np.arccos(dotProduct)/I3Units.deg)   

    alpha = [np.arccos(cosA)/I3Units.deg for cosA in cosAlpha]

    percent50Error = []
    percent90Error = []
    for errorList in listsError:
        errorList.sort()
        index50Per = int(0.5*len(errorList))
        index90Per = int(0.9*len(errorList))
    
        if len(errorList) == 0:
            percent50Error.append(-10)
            percent90Error.append(-10)
        else:
            percent50Error.append(errorList[index50Per])
            percent90Error.append(errorList[index90Per])

    return percent50Error, percent90Error, alpha, cosAlpha, unsuccessfulRecoEnergy, successfulRecoEnergy

binsE = np.linspace(3,7,11)
fitanalysisImproved = fitanalysis(infileImproved, binsE, "improved")
fitAnalysisLinefit = fitanalysis(infileLineFit, binsE, "linefit")

plt.hist(fitanalysisImproved[3], histtype = 'step', bins = 20)
plt.xlabel(r'$\cos{\alpha}$')
plt.title('Distribution of Relative Angle of Muon and its Reconstruction')

plt.figure()
plt.hist(fitanalysisImproved[3], histtype = 'step', bins = 20, log = True)
plt.xlabel(r'$\cos{\alpha}$')
plt.title('Distribution of Relative Angle of Muon and its Reconstruction')

plt.figure()
plt.hist(fitanalysisImproved[2], histtype = 'step', bins = 20)
plt.xlabel(r'$\alpha$')
plt.title('Distribution of Relative Angle of Muon and its Reconstruction')

plt.figure()
plt.hist(fitanalysisImproved[2], histtype = 'step', bins = 20, log = True)
plt.xlabel(r'$\alpha}$')
plt.title('Distribution of Relative Angle of Muon and its Reconstruction')

plt.figure()
plt.step(binsE[:-1], fitAnalysisLinefit[0], where = 'post', label = "50th percentile, linefit")
plt.step(binsE[:-1], fitanalysisImproved[0], where = 'post', label = "50th percentile, improved")
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel("Angular Difference (degrees)")
plt.title("Reconstuction Error - Successful Reco Only")
plt.legend()

plt.figure()
plt.step(binsE[:-1], fitAnalysisLinefit[1], where = 'post', label = "90th percentile, linefit")
plt.step(binsE[:-1], fitanalysisImproved[1], where = 'post', label = "90th percentile, improved")
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel("Angular Difference (degrees)")
plt.title("Reconstuction Error - Successful Reco Only")
plt.legend()

binsE = np.linspace(3,7,10)

h1, edges= np.histogram(fitanalysisImproved[4], bins = binsE)
h2, edges = np.histogram(fitanalysisImproved[5], bins = binsE)
for i in range(len(edges)-1):
    if h2[i] <=0.0000001:
        h2[i] = 1.0
        h1[i] = 0

plt.figure()
ratioHist = h1*1.0 / h2
plt.step(edges[:-1], ratioHist, where = 'post')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel("Fraction")
plt.title("Fraction of Failed Reconstructions")

plt.show()
