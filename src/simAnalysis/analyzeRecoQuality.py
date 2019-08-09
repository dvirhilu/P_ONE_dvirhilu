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
infilenoC = dataio.I3File('/home/dvir/workFolder/I3Files/linefitReco/HorizGeo/NuGen_linefitReco_HorizGeo_noCoincidence.i3.gz')
infileInitTime = dataio.I3File('/home/dvir/workFolder/I3Files/linefitReco/HorizGeo/NuGen_linefitReco_HorizGeo_testNewAlgorithm.i3.gz')

def findIndex(energy, bins):
    logE = np.log10(energy)

    for i in range(len(bins)):
        if logE < bins[i+1]:
            return i
    
    raise ValueError("energy not in range")


def fitanalysis(infile, binsE, fitType):
    cosAlpha = []
    listsError = [[] for i in range(10)]

    for frame in infile:
        primary = frame["NuGPrimary"]
        mctree = frame["I3MCTree"]
        muon = dataclasses.I3MCTree.first_child(mctree, primary)
        
        if fitType == "linefit":
            recoParticle = frame["LineFitRecoParticle"]
        else:
            recoParticle = frame["ImprovedRecoParticle"]

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

    return percent50Error, percent90Error, alpha, cosAlpha

binsE = np.linspace(3,7,11)
percent50ErrorImproved, percent90ErrorImproved, alphaImproved, cosAlphaImproved = fitanalysis(infile, binsE, "improved")
percent50ErrornoC, percent90ErrornoC, alphanoC, cosAlphanoC = fitanalysis(infilenoC, binsE, "linefit")
percent50ErrorC, percent90ErrorC, alphaC, cosAlphaC = fitanalysis(infileInitTime, binsE, "linefit")

'''
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
'''

print percent50ErrorC
print percent50ErrorImproved
plt.figure()
plt.step(binsE[:-1], percent50ErrorC, where = 'post', label = "50th percentile, with coincidence")
#plt.step(binsE[:-1], percent50ErrornoC, where = 'post', label = "50th percentile, no coincidence")
plt.step(binsE[:-1], percent50ErrorImproved, where = 'post', label = "50th percentile, improved")
plt.step(binsE[:-1], percent90ErrorC, where = 'post', label = "90th percentile, with coincidence")
#plt.step(binsE[:-1], percent90ErrornoC, where = 'post', label = "90th percentile, no coincidence")
plt.step(binsE[:-1], percent90ErrorImproved, where = 'post', label = "90th percentile, improved")
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel("Angular Difference (degrees)")
plt.title("LineFit Error")
plt.legend()

plt.figure()
plt.step(binsE[:-1], percent50ErrorC, where = 'post', label = "50th percentile, with coincidence")
#plt.step(binsE[:-1], percent50ErrornoC, where = 'post', label = "50th percentile, no coincidence")
plt.step(binsE[:-1], percent50ErrorImproved, where = 'post', label = "50th percentile, improved")
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel("Angular Difference (degrees)")
plt.title("LineFit Error")
plt.legend()

plt.figure()
plt.step(binsE[:-1], percent90ErrorC, where = 'post', label = "90th percentile, with coincidence")
#plt.step(binsE[:-1], percent90ErrornoC, where = 'post', label = "90th percentile, no coincidence")
plt.step(binsE[:-1], percent90ErrorImproved, where = 'post', label = "90th percentile, improved")
plt.xlabel(r'$log_{10}\, E/GeV$')
plt.ylabel("Angular Difference (degrees)")
plt.title("LineFit Error")
plt.legend()

plt.show()
