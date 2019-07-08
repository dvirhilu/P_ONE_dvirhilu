#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "Find neutrino distribution from NuGen simulation")
parser.add_argument('-i', '--infile', help = "input file" )
args = parser.parse_args()

# open file
infile = dataio.I3File(args.infile)

# get flux data
# make output file
infilePathStrings = args.infile.split("/")
infileName = infilePathStrings[len(infilePathStrings)-1]
infileAttributes = infileName.split("_")

fluxFileName = "fluxData_" + infileAttributes[0] + "_" + infileAttributes[2] + "_" + infileAttributes[3] + ".dat"
fluxFilePath = "/home/dvir/workFolder/P_ONE_dvirhilu/src/I3FileCode/simAnalysis/fluxData/"

fluxData = np.loadtxt(fluxFilePath + fluxFileName, unpack = True)
fluxDataMap = {eventID:fluxMult for eventID, fluxMult in zip(fluxData[0], fluxData[1])}

# get all Q frames
qframes = []
while infile.more():
    qframes.append(infile.pop_daq())

# get distribution of direction, position, and energy
zenith = []
azimuth = []
energy = []
x = []
y = []
z = []
weight = []
for frame in qframes:
    event_id = frame["I3EventHeader"].event_id
    primary = frame["NuGPrimary"]
    weightDict = frame["I3MCWeightDict"]
    oneWeight = weightDict["OneWeight"]
    numEvents = weightDict["NEvents"]
    fluxMult = fluxDataMap[event_id]
    weight.append(fluxMult*oneWeight/numEvents/2)

    zenith.append(np.cos(primary.dir.zenith))
    azimuth.append(np.cos(primary.dir.azimuth))
    energy.append(primary.energy)
    x.append(primary.pos.x)
    y.append(primary.pos.y)
    z.append(primary.pos.z)

logE = np.log10(energy)

plt.hist(logE, histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Log Energy Distribution")
plt.xlabel("LogE (Energy in GeV, log base 10)")
plt.ylabel("Number of Occurences")

plt.figure()
plt.hist(azimuth, histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Angular Distribution (Azimuth)")
plt.xlabel("Cosine of the Azimuth Angle")
plt.ylabel("Number of Occurences")

plt.figure()
plt.hist(zenith, histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Angular Distribution (Zenith)")
plt.xlabel("Cosine of the Zenith Angle")
plt.ylabel("Number of Occurences")

plt.figure()
plt.hist(x, histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Position Distribution (x)")
plt.xlabel("x Coordinate")
plt.ylabel("Number of Occurences")

plt.figure()
plt.hist(y, histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Position Distribution (y)")
plt.xlabel("y Coordinate")
plt.ylabel("Number of Occurences")

plt.figure()
plt.hist(z, histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Position Distribution (z)")
plt.xlabel("z Coordinate")
plt.ylabel("Number of Occurences")

plt.show()