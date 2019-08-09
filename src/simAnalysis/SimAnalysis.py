#!/usr/bin/env python

'''
SimAnalysis contains helper functions commonly used when writing
scripts to analyze results from simulation sets. Mainly, it allows
to easily decompose geometries into smaller sets by specifying a list
of omkeys, and allows to impose harder cuts on simulations so that different
geometries/cuts could be tried out without running the simulation again.

To import make sure that P_ONE_dvirhilu/src is in PYTHONPATH
'''

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame
import numpy as np
from numpy import linalg as la

# Calculates energy weight applied to frame based on predefined function
# 
# @Param:
# frame:            frame for the event in question
# energyWeighting:  the energy weighting function used (passed as a function
#                   and must is simply called once event info is parsed) 
# 
# @Return: 
# the energy weighting for the event
def weightE(frame, energyWeighting):
    primary = frame["NuGPrimary"]

    return energyWeighting(primary.energy)

# Calculates zenith weight applied to frame based on predefined function
# 
# @Param:
# frame:            frame for the event in question
# energyWeighting:  the zenith weighting function used (passed as a function
#                   and must is simply called once event info is parsed) 
# 
# @Return: 
# the zenith weighting for the event
def weightZenith(frame, zenithWeighting):
    primary = frame["NuGPrimary"]

    return zenithWeighting(primary.dir.zenith)

# Placeholder weighting function in case no weighting function is passed
# 
# @Param:
# weightingParam:   a dummy parameter which remains unsused. Since weighting 
#                   functions either take the particle energy or zenith for
#                   computation, a paramater is passed for compatability  
# 
# @Return: 
# 1 regardless of events (all events considered equal)
def defaultWeighter(weightingParam):
    return 1

# A function that determines whether a DOM passes the cut or not. The
# function checks whether a threshold number of hits was registered in a 20ns
# time window (any 20ns window counts)    
# 
# @Param:
# mcpeList:         The list of hits the DOM registered   
# hitThresh:        The number of hits required in the window for the DOM to
#                   be passed    
# 
# @Return: 
# A boolean variable indicating whether the DOM passed or not
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

# A function that determines whether a frame passes the cut or not. The
# function checks whether a threshold number of DOMs passed given a 
# threshold number of hits needed in the 20ns window
# 
# @Param:
# frame:            The frame in question   
# domsUsed:         List of omkeys used for the analysis. Allows to look 
#                   at smaller geometries from a larger geometry sim file
# hitThresh:        Hit threshold for a single DOM
# domThresh:        Number of passed DOMs needed to pass the frame 
# 
# @Return: 
# A boolean variable indicating whether the frame passed or not
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

# Calcualtes the number of retained frames given a list of unmasked DOMs,
# hit threshold to keep a DOM, and the number of DOMs to keep the frame.
# Two weighting functions, one for the primary's energy and one for the
# zenith angle can be provided to weight events based on desired particle
# characteristics to further customize results. This function is meant to
# serve as a quick first pass filter when testing different geometries from
# the dense one      
# 
# @Param:
# infileList:       List of I3Files containing simulation data   
# domsUsed:         List of omkeys used for the analysis. Allows to look 
#                   at smaller geometries from a larger geometry sim file
# hitThresh:        Hit threshold for a single DOM
# domThresh:        Number of passed DOMs needed to pass the frame
# energyWeighting:  A function (or functor) that returns the energy weighting
#                   given the particle's energy. Requires a callable object 
#                   that only a single parameter (particle's energy)  
# zenithWeighting:  A function (or functor) that returns the zenith weighting
#                   given the particle's zenith direction (keep in mind that 
#                   this refers to I3Particle.dir.zenith, which is the direction
#                   that the particle is coming from). Requires a callable object
#                    that only a single parameter (particle's zenith)  
# 
# @Return: 
# A number representing the weighted frame total
def calculateRetainedFrames(infileList, domsUsed, hitThresh, domThresh, energyWeighting = defaultWeighter, zenithWeighting = defaultWeighter):
    retainedFrames = 0
    for infile in infileList:
        infile.rewind()
        for frame in infile:
            if passFrame(frame, domsUsed, hitThresh, domThresh):
                retainedFrames += 1 * weightE(frame, energyWeighting) * weightZenith(frame, zenithWeighting)

    return retainedFrames

# A modified binary search algorithm to find all hits within a 20ns time 
# window. Takes entire list of MCPE hits, finds a 20ns window where the
# number of hits surpasses the hit threshold, and outputs those hits. If
# there are multiple windows, only the first occurence (in time) is taken. 
# NOTE: function assumes the list has already passed through the passDOM 
#       function. If the resulting list does not meet that criteria, the 
#       function raises a ValueError as it could cause issues for future 
#       functions relying on the significant MCPEs 
# 
# @Param:
# mcpeList:         The list of hits the DOM registered   
# hitThresh:        The number of hits required in the window for the DOM to
#                   be passed    
# 
# @Return: 
# An MCPESeries object with all the hits that were within the time window
def getSignificantMCPEs(mcpeList, hitThresh):
    
    if not passDOM(mcpeList, hitThresh):
        raise ValueError("The mcpeList given has not passed the passDOM function. Consider filtering with passDOM before running this function")
    
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
            # then search in timeList[middle + 1 to highest]
            if(timeList[middle] <= mcpe.time + 20): 
                lowest = middle + 1 
            else: 
            # else search in timeList[llowest to middle-1] 
                highest = middle - 1
        
        highIndex = highest

        lowest = 0
        while(lowest <= highest):
            middle = int((lowest+highest)/2)
            # if mcpe.time is greater than or equal to timeList[middle],  
            # then search in timeList[middle + 1 to highest]
            if(timeList[middle] < mcpe.time): 
                lowest = middle + 1 
            else: 
            # else search in timeList[lowest to middle-1] 
                highest = middle - 1
        
        lowIndex = highest + 1

        if highIndex - lowIndex >= hitThresh - 1:
            break
    
    for mcpe in mcpeList:
        if mcpe.time >= timeList[lowIndex] and mcpe.time <= timeList[highIndex]:
            significantMCPEList.append(mcpe)

    # sanity check
    if len(significantMCPEList) < hitThresh:
        raise ValueError("There's a bug in the code: mcpeList passed the passDOM function but resulted in a list too small")
    
    return significantMCPEList

# Converts the MCPESeriesMap to one that contains siginficant hits. Significant
# hits are defined to be ones where the DOM passed the passDOM function. The
# function then writes the new MCPESeriesMap to the frame before returning it
# NOTE: function assumes the list has already passed through the passFrame 
#       function. If the resulting map does not meet that criteria, the a 
#       ValueError is raised as it could cause issues for future functions 
#       relying on the significant MCPESeriesMap
#  
# @Param:
# frame:            The frame to be appended then returned   
# domsUsed:         The omkeys used for the analysis. Allows to look at
#                   smaller geometries from a larger geometry sim file
# hitThresh:        Hit threshold for a single DOM
# domThresh:        Number of DOMs passing the hit threshold to pass the
#                   frame  
# 
# @Return: 
# The frame inputted to the function after the significant MCPESeriesMap object
# was appended to it with the key word "MCPESerieMap_significant_hits" 
def writeSigHitsMapToFrame(frame, domsUsed, hitThresh, domThresh):
    if not passFrame(frame, domsUsed, hitThresh, domThresh):
        raise ValueError("The mcpeList given has not passed the passFrame function. Consider filtering with passFrame before running this function")
    
    mcpeMap = frame["MCPESeriesMap"]
    significantMCPEMap = simclasses.I3MCPESeriesMap()
    for omkey in domsUsed:
        if omkey not in mcpeMap:
            continue
        if passDOM(mcpeMap[omkey], hitThresh):
            significantMCPEMap[omkey] = getSignificantMCPEs(mcpeMap[omkey], hitThresh)  
    
    # sanity check
    if len(mcpeMap) < domThresh:
        raise ValueError("There's a bug in the code: mcpeMap passed the passFrame function but resulted in a map too small")
    
    frame.Put("MCPESeriesMap_significant_hits", significantMCPEMap)

    return frame

# parses the hit information requires for the linefit reconstruction and returns
# a list of tuples containing the required hit information
# NOTE: function assumes the frame already contains the significant MCPESeriesMap
#       and tries to call on the "MCPESeriesMap_significant_hits" key. If the key
#       does not exist, the function will raise a ValueError. It is much more  
#       efficient to have the significant hits MCPESeriesMap to use rather than
#       having to search for the significant hits every time, so most methods rely
#       on the frame key being there.   
#
# @Param:
# frame:            The frame containing the event information   
# geometry:         An I3Geometry object from the gcd file used to produce the
#                   simulation data 
# hitThresh:        Hit threshold for a single DOM
#
# @Return: 
# A list of tuples. Each tuple contains information for a single hit (x,y,z,t)
def getLinefitDataPoints(frame, geometry):
    if not frame.Has("MCPESeriesMap_significant_hits"):
        raise ValueError("Frame does not contain the MCPESeriesMap_significant_hits. Make sure to run the writeSigHitsToFrame function")

    mcpeMap = frame["MCPESeriesMap_significant_hits"]
    geoMap = geometry.omgeo
    data = []
    for omkey, mcpeList in mcpeMap:
        timeList = [mcpe.time for mcpe in mcpeList]
        time = min(timeList)
        position = geoMap[omkey].position
        for i in range(len(timeList)):
            data.append( (position.x, position.y, position.z, time) )

    return data

# Given a list of x and y points, computes the parameters of a least squares fit
# line
#  
# @Param:
# x:                The x component of points in the fit   
# y:                The y component of points in the fit
# @Return: 
# A tuple containing the slope and y intercept of least squares fit line
def fitLeastSquaresLine(x, y):
    # for a given system X*c = y 
    # where X is a truncated Vandermond matrix (truncated in this case to be degree 1)
    #       c is the vector containing the polyniomial coefficients
    #       y is the vector containing all y values
    # The least squares solution is given as c = (X.T * X)^(-1) * (X.T * y)
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

# Given datapoints containing the hit informations (from getLineFitDatapoints),
# the function computes the parameters for the linefit reconstruction of the 
# particle 
# 
# @Param:
# datapoints:       A list of tuples. Each tuple contains hit information in the
#                   format (x,y,z,t) 
# 
# @Return: 
# A tuple containing the reconstructed particle's information. This is in the
# of an I3Direction object for the particle's direction, a double for the particle's
# speed, and an I3Position object for the particle vertex (position at t=0)  
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