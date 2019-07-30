#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units, I3Frame
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib

def weightE(frame, energyWeighting):
    primary = frame["NuGPrimary"]

    return energyWeighting(primary.energy)

def weightZenith(frame, zenithWeighting):
    primary = frame["NuGPrimary"]

    return zenithWeighting(primary.zenith)

def defaultWeighter(weightingParameter):
    return 1


def passDOM(mcpeList, hitThresh):
    if len(mcpeList) < hitThresh:
        return False

    timeList = [mcpe.time for mcpe in mcpeList]

    if timeList[len(timeList) - 1] - timeList[0] < 20:
        return True

    for i in range(len(timeList) - hitThresh) 
    
def passFrame(frame, domsUsed, hitThresh, domThresh):
    if frame.Stop != I3Frame.DAQ:
        return False
    mcpeMap = frame["MCPESeriesMap"]
    
    domCount = 0
    for dom in domsUsed:
        if passDOM(mcpeMap[dom], hitThresh):
            domCount += 1
    
    if domCount < domThresh:
        return False
    
    return True

def calculateRetainedFrames(infileList, domsUsed, hitThresh, domThresh, energyWeighting = defaultWeighter, zenithWeighting = defaultWeighter):
    totalFrames = 0
    retainedFrames = 0
    frameList = []
    for infile in infileList:
        for frame in infile:
            totalFrames += 1
            if passFrame(frame, domsUsed, hitThresh, domThresh):
                retainedFrames += 1 * weightE(frame, energyWeighting) * weightZenith(frame, zenithWeighting)
                frameList.append(frame)
    

    return frameList, retainedFrames, 
