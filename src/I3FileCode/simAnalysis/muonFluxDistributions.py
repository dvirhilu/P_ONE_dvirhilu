#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units
import matplotlib.pyplot as plt
import numpy as np

# open file
infile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/muongun/muongun_step1/MuonGun_step1_139005_000000.i3.bz2')
outfile = open('muonDirectionz.txt', 'w+')

# get all Q frames
qframes = []
while infile.more():
    qframes.append(infile.pop_daq())

# get distribution of direction and energy
zenith = []
azimuth = []
energy = []
zdirection = []
i=1
for frame in qframes:
    # get propagated muon
    mctree = frame["I3MCTree"]
    primary = mctree.primaries
    muon = dataclasses.I3MCTree.first_child(mctree, primary[0].id)

    # add muon energy
    energy.append(muon.energy / I3Units.GeV)

    # get muon direction
    muonDirection = muon.dir
    azimuth.append(np.cos(muonDirection.azimuth))
    zenith.append(np.cos(muonDirection.zenith))
    outfile.write("direction for muon number %d: \n" % (i))
    outfile.write("x: " + str(muonDirection.x) + "\n")
    outfile.write("y: " + str(muonDirection.y) + "\n")
    outfile.write("z: " + str(muonDirection.z) + "\n")
    outfile.write("azimuth: " + str(muonDirection.azimuth / I3Units.deg) + "\n")
    outfile.write("zenith: " + str(muonDirection.zenith / I3Units.deg) + "\n")
    outfile.write("phi: " + str(muonDirection.phi / I3Units.deg) + "\n")
    outfile.write("theta: " + str(muonDirection.theta / I3Units.deg) + "\n\n\n")
    i += 1


# plot histograms
xmin = 500 * I3Units.GeV
xmax = 1 * I3Units.TeV
plt.hist(energy, range = (xmin,xmax), histtype = "step", log = True, bins = 25)
plt.title("Muon Energy Distribution")
plt.xlabel("E(GeV)")
plt.show()

plt.figure()
xmin = -1
xmax = 1
plt.hist(azimuth, range = (xmin,xmax), histtype = "step", log = True, bins = 45)
plt.title("Muon Angular Distribution (Azimuth)")
plt.xlabel("Azimuth (degrees)")
plt.show()

plt.figure()
xmin = -1
xmax = 1
plt.hist(zenith, range = (xmin,xmax), histtype = "step", log = True, bins = 45)
plt.title("Muon Angular Distribution (Zenith)")
plt.xlabel("Zenith (degrees)")
plt.show()

outfile.close()