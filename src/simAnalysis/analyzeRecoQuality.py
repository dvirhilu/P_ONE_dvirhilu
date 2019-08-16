#!/usr/bin/env python

from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame
from icecube.dataclasses import I3Particle
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from simAnalysis.SimAnalysis import makeRatioHist

infile5LineGeoLineFit = dataio.I3File('/home/dvir/workFolder/I3Files/linefitReco/partialDenseGeo/NuGen_linefitReco_partialDenseGeo_5LineGeo.i3.gz')
infile10LineGeoLineFit = dataio.I3File('/home/dvir/workFolder/I3Files/linefitReco/partialDenseGeo/NuGen_linefitReco_partialDenseGeo_10LineGeo.i3.gz')
infile5LineGeoImproved = dataio.I3File('/home/dvir/workFolder/I3Files/improvedReco/partialDenseGeo/NuGen_improvedReco_partialDenseGeo_5LineGeo.i3.gz')
infile10LineGeoImproved = dataio.I3File('/home/dvir/workFolder/I3Files/improvedReco/partialDenseGeo/NuGen_improvedReco_partialDenseGeo_10LineGeo.i3.gz')

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
fitanalysis5LineLF = fitanalysis(infile5LineGeoLineFit, binsE, "linefit")
fitAnalysis10LineLF = fitanalysis(infile10LineGeoLineFit, binsE, "linefit")
fitanalysis5LineIF = fitanalysis(infile5LineGeoImproved, binsE, "improved")
fitAnalysis10LineIF = fitanalysis(infile10LineGeoImproved, binsE, "improved")

plt.hist(fitanalysis5LineLF[2], log = True, histtype = 'step', bins = 20, label = '5 line, linefit')
plt.hist(fitAnalysis10LineLF[2], log = True, histtype = 'step', bins = 20, label = '10 line, linefit')
plt.hist(fitanalysis5LineIF[2], log = True, histtype = 'step', bins = 20, label = '5 line, chi-squared')
plt.hist(fitAnalysis10LineIF[2], log = True, histtype = 'step', bins = 20, label = '10 line, chi-squared')
plt.xlabel(r'$\alpha$')
plt.title('Distribution of Relative Angle of Muon and its Reconstruction')
plt.legend()

plt.figure()
plt.step(binsE[:-1], fitanalysis5LineLF[0], where = 'post', label = "5line, linefit")
plt.step(binsE[:-1], fitAnalysis10LineLF[0], where = 'post', label = "10line, linefit")
plt.step(binsE[:-1], fitanalysis5LineIF[0], where = 'post', label = " 5line, chi-squared")
plt.step(binsE[:-1], fitAnalysis10LineIF[0], where = 'post', label = "10line, chi-squared")
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel("Angular Difference (degrees)")
plt.title("Reconstuction Error - Successful Reco Only, 50th Percentile")
plt.legend()

plt.figure()
plt.step(binsE[:-1], fitanalysis5LineLF[1], where = 'post', label = "5line, linefit")
plt.step(binsE[:-1], fitAnalysis10LineLF[1], where = 'post', label = "10line, linefit")
plt.step(binsE[:-1], fitanalysis5LineIF[1], where = 'post', label = " 5line, chi-squared")
plt.step(binsE[:-1], fitAnalysis10LineIF[1], where = 'post', label = "10line, chi-squared")
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel("Angular Difference (degrees)")
plt.title("Reconstuction Error - Successful Reco Only, 90th Percentile")
plt.legend()

'''
binsE = np.linspace(3,7,10)

plt.figure()
ratioHist, edges = makeRatioHist(fitanalysisImproved[4], fitanalysisImproved[5], weights1 = np.ones(len(fitanalysisImproved[4])), weights2 = np.ones(len(fitanalysisImproved[5])), bins = binsE)
plt.step(edges[:-1], ratioHist, where = 'post')
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel("Fraction")
plt.title("Fraction of Failed Reconstructions")
'''
plt.show()
