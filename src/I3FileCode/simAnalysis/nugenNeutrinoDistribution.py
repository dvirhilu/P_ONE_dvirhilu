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
ptypeMap = {}
zenith500prim = []
zenith500sec = []
weight500 = []
zenithZ300prim = []
zenithZ300sec = []
weightZ300 = []
angDiff = []

for frame in qframes:
    event_id = frame["I3EventHeader"].event_id
    primary = frame["NuGPrimary"]
    eventWeight = frame["EventWeight"].value
    secondaryParticle = dataclasses.I3MCTree.children(frame["I3MCTree"], primary)[0]
    weight.append(eventWeight)
    zenith.append(np.cos(secondaryParticle.dir.zenith))
    azimuth.append(np.cos(secondaryParticle.dir.azimuth))
    energy.append(primary.energy)
    x.append(secondaryParticle.pos.x)
    y.append(secondaryParticle.pos.y)
    z.append(secondaryParticle.pos.z)

    ptype = secondaryParticle.type
    if ptype in ptypeMap:
        ptypeMap[ptype] += 1
    else:
        ptypeMap[ptype] = 1

    if secondaryParticle.pos.z > 300:
        zenithZ300prim.append(np.cos(primary.dir.zenith))
        zenithZ300sec.append(np.cos(secondaryParticle.dir.zenith))
        weightZ300.append(weight[len(weight)-1])

    if x[len(x)-1]**2 + y[len(y)-1]**2 +z[len(z)-1]**2 < 500**2:
        zenith500sec.append(np.cos(secondaryParticle.dir.zenith))
        weight500.append(weight[len(weight) - 1])
        zenith500prim.append(np.cos(primary.dir.zenith))

    angDiff.append(primary.dir.x*secondaryParticle.dir.x + primary.dir.y*secondaryParticle.dir.y + primary.dir.z*secondaryParticle.dir.z)

logE = np.log10(energy)
r = [np.sqrt(x[i]**2 + y[i]**2 + z[i]**2) for i in range(len(x))]
farEvents = [value for value in r if value > 5000]

plt.hist(logE, histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Log Energy Distribution")
plt.xlabel("LogE (Energy in GeV, log base 10)")

plt.figure()
plt.hist(azimuth, histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Angular Distribution (Azimuth)")
plt.xlabel("Cosine of the Azimuth Angle")

plt.figure()
plt.hist(zenith, histtype = "step", log = True, weights = weight, bins = 20)
plt.title("Weighted Muon Angular Distribution (Zenith)")
plt.xlabel("Cosine of the Zenith Angle")

plt.figure()
plt.hist(x, histtype = "step", log = True, weights = weight, bins = 50)
plt.title("Weighted Muon Position Distribution (x)")
plt.xlabel("x Coordinate")

plt.figure()
plt.hist(y, histtype = "step", log = True, weights = weight, bins = 50)
plt.title("Weighted Muon Position Distribution (y)")
plt.xlabel("y Coordinate")

plt.figure()
plt.hist(z, histtype = "step", log = True, weights = weight, bins = 50)
plt.title("Weighted Muon Position Distribution (z)")
plt.xlabel("z Coordinate")

plt.figure()
plt.hist(r, histtype = "step", log = True, weights = weight, bins = 50)
plt.title("Weighted Muon Distance Distribution")
plt.xlabel("Distance from 0")
'''
plt.figure()
plt.hist(zenith500prim, histtype = "step", log = True, weights = weight500, bins = 20)
plt.title("Weighted Neutrino Angular Distribution (Zenith), Resulting Muon < 500m Away")
plt.xlabel("Cosine of the Zenith Angle")

plt.figure()
plt.hist(zenith500sec, histtype = "step", log = True, weights = weight500, bins = 20)
plt.title("Weighted Muon Angular Distribution (Zenith), < 500m Away")
plt.xlabel("Cosine of the Zenith Angle")

plt.figure()
plt.hist(zenithZ300prim, histtype = "step", log = True, weights = weightZ300, bins = 20)
plt.title("Weighted Neutrino Angular Distribution (Zenith), Resulting Muon Z>300m")
plt.xlabel("Cosine of the Zenith Angle")

plt.figure()
plt.hist(zenithZ300sec, histtype = "step", log = True, weights = weightZ300, bins = 20)
plt.title("Weighted Muon Angular Distribution (Zenith), Z>300m")
plt.xlabel("Cosine of the Zenith Angle")

plt.figure()
plt.hist(angDiff, histtype = "step", log = True, bins = 100)
plt.title("Weighted Angular Difference Distribution Neutrino vs. Muon")
plt.xlabel("Cosine of the Angular Difference")
'''
print(len(zenith500sec))
print(len(farEvents))
print(ptypeMap)
plt.show()