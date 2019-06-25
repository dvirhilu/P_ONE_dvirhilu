#!/usr/bin/env python

from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey
import matplotlib.pyplot as plt
import numpy as np

infile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/muongun/muongun_step3/MuonGun_step3_139005_000000.i3.bz2')

# get all Q frames
qframes = []
while(infile.more()):
    qframes.append(infile.pop_daq())

# create a list of number of photons detected
eventNum = np.linspace(1, len(qframes), len(qframes))
peCount = np.array([])
domHits = np.array([])
nDomsWLight = np.array([])

# fill information about hits
for frame in qframes:
    series_map = frame["MCPESeriesMap"]
    nDomsWLight = np.append(nDomsWLight, len(series_map) )
    hits = 0
    pe = 0
    
    for omkey in series_map.keys():
        hitArray = series_map[omkey]
        hits += len(hitArray)
        for hit in hitArray:
            pe += hit.npe
    
    domHits = np.append(domHits, hits)
    peCount = np.append(peCount, pe)

    


# plot hit information vs. event number
plt.hist(peCount, color = 'blue', label = 'total PE', histtype = 'step', bins = 50)
plt.hist(domHits, color = 'green', label = 'total DOM hits', histtype = 'step', bins = 50)
plt.hist(nDomsWLight, color = 'red', label = 'number of DOMs that were hit', histtype = 'step', bins = 50)
plt.title("Hit Information for Each Event")
plt.xlabel("Event Number")
plt.ylabel("Number of Occurences")
plt.legend(loc = 'upper right')
plt.show()