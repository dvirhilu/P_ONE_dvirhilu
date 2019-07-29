#!/usr/bin/env python

from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame
from icecube.dataclasses import I3Particle
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as la
import argparse

parser = argparse.ArgumentParser(description = "Creates a reconstruction of the muon track using a linear least squares fit on the pulses")
parser.add_argument( '-n', '--minFileNum', help = "smallest file number used" )
parser.add_argument( '-N', '--maxFileNum', help = "largest file number used")
parser.add_argument( '-H', '--hitThresh', help = "threshold of hits for the DOM to be considered")
parser.add_argument( '-d', '--domThresh', help = "threshold of hit DOMs for the frame to be considered")
parser.add_argument( '-g', '--GCDType', help = "type of geometry used for the simulation set")
args = parser.parse_args()

if args.GCDType == 'testString':
    gcdPath = '/home/dvir/workFolder/I3Files/gcd/testStrings/HorizTestString_n15_b100.0_v50.0_l1_simple_spacing.i3.gz'
elif args.GCDType == 'HorizGeo':
    gcdPath = '/home/dvir/workFolder/I3Files/gcd/corHorizgeo/CorrHorizGeo_n15_b100.0_a18.0_l3_rise_fall_offset_simple_spacing.i3.gz'
elif args.GCDType == 'IceCube':
    gcdPath = '/home/dvir/workFolder/I3Files/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz'
elif args.GCDType == 'cube':
    gcdPath = '/home/dvir/workFolder/I3Files/gcd/cube/cubeGeometry_1600_15_50.i3.gz'
else:
    raise RuntimeError("Invalid GCD Type")

infileList = []
for i in range(int(args.minFileNum), int(args.maxFileNum)+1):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/' + str(args.GCDType) + '/NuGen_step3_' + str(args.GCDType) + '_' + str(i) + '.i3.gz')
    infileList.append(infile)

outfile = dataio.I3File('/home/dvir/workFolder/I3Files/linefitReco/'+ str(args.GCDType) + '/NuGen_linefitReco_' + str(args.GCDType) + '_noCoincidence.i3.gz', 'w')
gcdfile = dataio.I3File(gcdPath)
geometry = gcdfile.pop_frame()["I3Geometry"]

def passFrame(frame, domThresh, hitThresh):
    if frame.Stop != I3Frame.DAQ:
        return False
    mcpeMap = frame["MCPESeriesMap"]
    domCount = 0
    for mcpeList in mcpeMap.values():
        if len(mcpeList) >= hitThresh:
            domCount += 1
    
    if domCount < domThresh:
        return False
    
    return True

def getRecoDataPoints(frame, geometry, hitThresh):
    mcpeMap = frame["MCPESeriesMap"]
    geoMap = geometry.omgeo
    significantMCPEMap = simclasses.I3MCPESeriesMap()
    data = []
    for omkey, mcpeList in mcpeMap:
        if len(mcpeList) >= hitThresh:
            significantMCPEMap[omkey] = mcpeList
            for mcpe in mcpeList:
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

for infile in infileList:
    for frame in infile:
        if passFrame(frame, int(args.domThresh), int(args.hitThresh)):
            datapoints = getRecoDataPoints(frame, geometry, int(args.hitThresh))
            direction, speed, vertex = reconstructParticleParams(datapoints)

            recoParticle = I3Particle()
            recoParticle.shape = I3Particle.InfiniteTrack
            recoParticle.fit_status = I3Particle.OK
            recoParticle.dir = direction
            recoParticle.speed = speed
            recoParticle.pos = vertex
            recoParticle.time = 0

            frame.Put("LineFitRecoParticle", recoParticle)
            outfile.push(frame)


outfile.close()