#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray, NuFlux
from icecube
from icecube.icetray import I3Units
import matplotlib.pyplot as plt
import numpy as np

# open file
infile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/nugen/nugenStep1/NuGen_step1_testString_000899_.i3.gz')

# flux model
flux = NuFlux.makeFlux('honda2006').getFlux

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
    primary = frame["NuGPrimary"]
    weightDict = frame["I3MCWeightDict"]
    oneWeight = weightDict["OneWeight"]
    numEvents = weightDict["NEvents"]
    fluxMult = flux(primary.type, primary.energy, np.cos(primary.dir.zenith))
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