#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as la
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
    timeList.sort()

    if timeList[len(timeList) - 1] - timeList[0] < 20:
        return True

    for i in range(len(timeList) - hitThresh):
        if timeList[i+hitThresh-1] - timeList[i] < 20:
            return True 
        
    return False

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
            totalFrames += 1 * weightE(frame, energyWeighting) * weightZenith(frame, zenithWeighting)
            if passFrame(frame, domsUsed, hitThresh, domThresh):
                retainedFrames += 1 * weightE(frame, energyWeighting) * weightZenith(frame, zenithWeighting)
                frameList.append(frame)
    

    return frameList, retainedFrames, totalFrames


def getSignificantMCPEs(mcpeList, hitThresh):
    for mcpe in mcpeList:


def getRecoDataPoints(frame, geometry, hitThresh):
    mcpeMap = frame["MCPESeriesMap"]
    geoMap = geometry.omgeo
    significantMCPEMap = simclasses.I3MCPESeriesMap()
    data = []
    for omkey, mcpeList in mcpeMap:
        if passDOM(omkey, hitThresh):
            significantMCPEMap[omkey] = getSignificantMCPEs(mcpeList, hitThresh)
            for mcpe in significantMCPEMap[omkey]:
                time = mcpe.time
                position = geoMap[omkey].position
                data.append([position.x, position.y, position.z, time])
    
    frame.Put("MCPESeriesMap_significant_hits", significantMCPEMap)

    return data

def fitLeastSquaresLine(x, y):
    xMatrix = np.column_stack( (x, np.ones(len(x))) )
    yVector = np.array(y).T
    leastSquaresMatrix = np.matmul(xMatrix.T, xMatrix)
    leastSquaresVector = np.matmul(xMatrix.T, yVector)
    fitCoefficients = np.matmul( la.inv(leastSquaresMatrix), leastSquaresVector)

    slope = fitCoefficients[0]
    intercept = fitCoefficients[1]

    return slope, intercept

def reconstructParticleParams(datapoints):
    x = [data[0] for data in datapoints]
    y = [data[1] for data in datapoints]
    z = [data[2] for data in datapoints]
    t = [data[3] for data in datapoints]

    xVelocity, x = fitLeastSquaresLine(t, x)
    yVelocity, y = fitLeastSquaresLine(t, y)
    zVelocity, z = fitLeastSquaresLine(t, z)

    direction = dataclasses.I3Direction(xVelocity, yVelocity, zVelocity)
    speed = np.sqrt(xVelocity**2 + yVelocity**2 + zVelocity**2)
    vertex = dataclasses.I3Position(x,y,z)

    return direction, speed, vertex