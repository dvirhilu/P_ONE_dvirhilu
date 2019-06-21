#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units
import matplotlib.pyplot as plt

# open file
infile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/muongun/muongun_step1/MuonGun_step1_139005_000000.i3.bz2')

# get all Q frames
qframes = []
while infile.more():
    qframes.append(infile.pop_daq())

# get distribution of direction and energy
zenith = []
azimuth = []
energy = []
for frame in qframes:
    # get propagated muon
    mctree = frame["I3MCTree"]
    primary = mctree.primaries
    muon = dataclasses.I3MCTree.first_child(mctree, primary[0].id)

    # add muon energy
    energy.append(muon.energy / I3Units.GeV)

    # get muon direction
    muonDirection = muon.dir
    azimuth.append(muonDirection.azimuth)
    zenith.append(muonDirection.zenith)


# plot histograms
xmin = 500 * I3Units.GeV
xmax = 1 * I3Units.TeV
plt.hist(energy, range = (xmin,xmax), histtype = "step", log = True, bins = 50)
plt.title("Muon Energy Distribution")
plt.xlabel("E(GeV)")
plt.show()

plt.figure()
xmin = 0
xmax = 360
plt.hist(azimuth, range = (xmin,xmax), histtype = "step", log = True, bins = 100)
plt.title("Muon Angular Distribution (Azimuth)")
plt.xlabel("Azimuth (degrees)")
plt.show()

plt.figure()
xmin = 0
xmax = 10
plt.hist(zenith, range = (xmin,xmax), histtype = "step", log = True, bins = 100)
plt.title("Muon Angular Distribution (Zenith)")
plt.xlabel("Zenith (degrees)")
plt.show()