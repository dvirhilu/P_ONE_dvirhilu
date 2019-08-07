#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame
import numpy as np
from numpy import linalg as la

def weightE(frame, energyWeighting):
    primary = frame["NuGPrimary"]

    return energyWeighting(primary.energy)

def weightZenith(frame, zenithWeighting):
    primary = frame["NuGPrimary"]

    return zenithWeighting(primary.dir.zenith)

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
        if dom not in mcpeMap:
            continue
        if passDOM(mcpeMap[dom], hitThresh):
            domCount += 1
    
    if domCount < domThresh:
        return False
    
    return True

def calculateRetainedFrames(infileList, domsUsed, hitThresh, domThresh, energyWeighting = defaultWeighter, zenithWeighting = defaultWeighter):
    totalFrames = 0
    retainedFrames = 0
    for infile in infileList:
        infile.rewind()
        for frame in infile:
            totalFrames += 1 * weightE(frame, energyWeighting) * weightZenith(frame, zenithWeighting)
            if passFrame(frame, domsUsed, hitThresh, domThresh):
                retainedFrames += 1 * weightE(frame, energyWeighting) * weightZenith(frame, zenithWeighting)

    return retainedFrames, totalFrames


def getSignificantMCPEs(mcpeList, hitThresh):
    timeList = [mcpe.time for mcpe in mcpeList]
    timeList.sort()
    highIndex = 0
    lowIndex = 0
    significantMCPEList = simclasses.I3MCPESeries()

    for mcpe in mcpeList:

        lowest = 0
        highest = len(timeList)-1
        while(lowest <= highest):
            middle = int((lowest+highest)/2)
            # if mcpe.time is greater than or equal to timeList[middle],  
            # then search in timeList[mid + 1 to h]
            if(timeList[middle] <= mcpe.time + 20): 
                lowest = middle + 1 
            else: 
            # else search in timeList[l to mid-1] 
                highest = middle - 1
        
        highIndex = highest

        lowest = 0
        while(lowest <= highest):
            middle = int((lowest+highest)/2)
            # if mcpe.time is greater than or equal to timeList[middle],  
            # then search in timeList[mid + 1 to h]
            if(timeList[middle] < mcpe.time): 
                lowest = middle + 1 
            else: 
            # else search in timeList[l to mid-1] 
                highest = middle - 1
        
        lowIndex = highest + 1

        if highIndex - lowIndex >= hitThresh - 1:
            break
    
    for mcpe in mcpeList:
        if mcpe.time >= timeList[lowIndex] and mcpe.time <= timeList[highIndex]:
            significantMCPEList.append(mcpe)
    
    return significantMCPEList

def getRecoDataPoints(frame, geometry, hitThresh):
    mcpeMap = frame["MCPESeriesMap"]
    geoMap = geometry.omgeo
    significantMCPEMap = simclasses.I3MCPESeriesMap()
    data = []
    for omkey, mcpeList in mcpeMap:
        if passDOM(mcpeList, hitThresh):
            significantMCPEMap[omkey] = getSignificantMCPEs(mcpeList, hitThresh)
            timeList = [mcpe.time for mcpe in significantMCPEMap[omkey]]
            time = min(timeList)
            position = geoMap[omkey].position
            for i in range(len(timeList)):
                data.append([position.x, position.y, position.z, time])
    
    frame.Put("MCPESeriesMap_significant_hits", significantMCPEMap)

    return data

def fitLeastSquaresLine(x, y):
    xMatrix = np.column_stack( (x, np.ones(len(x))) )
    yVector = np.array(y).T
    leastSquaresMatrix = np.matmul(xMatrix.T, xMatrix)
    if la.det(leastSquaresMatrix) == 0:
        print leastSquaresMatrix, xMatrix
    leastSquaresVector = np.matmul(xMatrix.T, yVector)
    fitCoefficients = np.matmul( la.inv(leastSquaresMatrix), leastSquaresVector)

    slope = fitCoefficients[0]
    intercept = fitCoefficients[1]

    return slope, intercept

def linefitParticleParams(datapoints):
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