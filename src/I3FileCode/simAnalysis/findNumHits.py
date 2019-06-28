#!/usr/bin/env python

from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey
import matplotlib.pyplot as plt
import numpy as np

step3infile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/muongun/muongun_step3/MuonGun_step3_139005_000000.i3.bz2')
custominfile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/muongun/customGenHitsMuongun/MuonGun_customGenHits_0.i3.gz')
missingHitsFile = open("MissingHits.txt", 'w')

# get all Q frames
step3qframes = []
while(step3infile.more()):
    step3qframes.append(step3infile.pop_daq())

customqframes = []
while(custominfile.more()):
    customqframes.append(custominfile.pop_daq())

# create a dictionary of event header and frame to ensure events match up
evtFrameMap = {}
for frame in customqframes:
    header = frame["I3EventHeader"]
    evtFrameMap[header.event_id] = frame

def getNumHitsInFrame(frame):
    mcpeSeriesMap = frame["MCPESeriesMap"]
    numHits = 0

    for omkey in mcpeSeriesMap.keys():
        hits = mcpeSeriesMap[omkey]
        numHits += len(hits)
    
    return numHits

step3Hits = []
customHits = []
event_id = []

for frame in step3qframes:
    header = frame["I3EventHeader"]
    evt_id = header.event_id
    if evt_id in evtFrameMap.keys():
        event_id.append(evt_id)
        customFrame = evtFrameMap[evt_id]
        step3Hits.append(getNumHitsInFrame(frame))
        customHits.append(getNumHitsInFrame(customFrame))
    else:
        hits = getNumHitsInFrame(frame)
        missingHitsFile.write("Event only in step3 file: \nEventID: " + str(evt_id) + ", Number of Hits: " + str(hits) + "\n\n")


for evt_id in evtFrameMap.keys():
    if not evt_id in event_id:
        hits = getNumHitsInFrame(evtFrameMap[evt_id])
        missingHitsFile.write("Event only in custom file: \nEventID: " + str(evt_id) + ", Number of Hits: " + str(hits) + "\n\n" )

# plot hit information vs. event number
plt.plot(event_id, step3Hits, 'r', label = "Number of Hits in Step 3")
plt.plot(event_id, customHits, 'b', label = "Number of Hits in Custom HitGen")
plt.title("Number of Hits in Each Event")
plt.xlabel("Event ID")
plt.ylabel("Total Number of Hits")
plt.legend()
plt.show()