#!/usr/bin/env python

from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, I3Frame
from icecube.dataclasses import I3Particle
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as la
import argparse, SimAnalysis

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

outfile = dataio.I3File('/home/dvir/workFolder/I3Files/linefitReco/'+ str(args.GCDType) + '/NuGen_linefitReco_' + str(args.GCDType) + '_testAlgorithm.i3.gz', 'w')
gcdfile = dataio.I3File(gcdPath)
geometry = gcdfile.pop_frame()["I3Geometry"]

for infile in infileList:
    for frame in infile:
        if SimAnalysis.passFrame(frame, geometry.omgeo.keys(), int(args.domThresh), int(args.hitThresh)):
            datapoints = SimAnalysis.getRecoDataPoints(frame, geometry, int(args.hitThresh))
            direction, speed, vertex = SimAnalysis.reconstructParticleParams(datapoints)

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