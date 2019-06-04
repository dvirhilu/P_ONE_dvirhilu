#!/usr/bin/env python

from icecube import dataio, dataclasses
import numpy as np
import matplotlib.pyplot as plt
import math

# initialize files
geofile = dataio.I3File('/home/dvir/workFolder/I3Files/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz')
infile = dataio.I3File('/home/dvir/workFolder/I3Files/Level2_IC86.2011_corsika.010281.001664.00.i3.bz2')
outfileC = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/generated/I3CascadeFile.I3.gz','w')
outfileM = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/generated/I3MuonFile.I3.gz', 'w')

# filter for cascades and output to oufileC
for frame in infile:
    fmask = frame["FilterMask"]
    cFilter = fmask["CascadeFilter_11"]
    if cFilter.condition_passed:
        outfileC.push(frame)

outfileC.close()

# rewind frame stack
infile.rewind()

# filter for muons and output to outfileM
for frame in infile:
    fmask = frame["FilterMask"]
    mFilter = fmask["MuonFilter_11"]
    if mFilter.condition_passed:
        outfileM.push(frame)
outfileM.close()

# rewind frame stack
infile.rewind()

# generate list of OfflinePulses and total charge for every event in the file
totPulses = []
totCharge = []
for frame in infile:
    offPulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, "OfflinePulses")
    pulseList = [pulse for omkey,pulses in offPulses for pulse in pulses]
    chargeList = [pulse.charge for pulse in pulseList]
    totCharge.append(sum(chargeList))
    totPulses.append(len(offPulses))

plt.hist(totPulses, histtype = "step", bins = 1000)
plt.xlabel("Number of DOMs That Saw a Photon in Each Event")
plt.ylabel("Number of Occurences")
plt.show()
plt.figure()
plt.hist(totCharge, histtype = "step", bins = 1000)
plt.xlabel("Total Charge in Each Event")
plt.ylabel("Number of Occurences")
plt.show()
