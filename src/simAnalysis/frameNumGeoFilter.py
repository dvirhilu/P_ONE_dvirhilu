#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame
import numpy as np
from simAnalysis import SimAnalysis
from experimentModelCode import FunctionClasses
import argparse
from iminuit import Minuit

'''
A script that looks at the number of frames retained by different decompositions of
a dense geometry. This allows to quickly discard geometries that are not of sufficient
quality without having to run a full simulation set for them.
'''

parser = argparse.ArgumentParser(description = "Creates a reconstruction of the muon track using a linear least squares fit on the pulses")
parser.add_argument( '-n', '--minFileNum', help = "smallest file number used" )
parser.add_argument( '-N', '--maxFileNum', help = "largest file number used")
parser.add_argument( '-H', '--hitThresh', help = "threshold of hits for the DOM to be considered")
parser.add_argument( '-D', '--domThresh', help = "threshold of hit DOMs for the frame to be considered")
parser.add_argument( '-g', '--GCDType', help = "type of geometry used for the simulation set")
parser.add_argument( '-l', '--numLayers', help = "number of layers the decomposed geometry will have")
parser.add_argument( '-o', '--outputName', help = "output name chosen for this specific iteration of the geometry")
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
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/' + str(args.GCDType) + '/NuGen_step3_' + 
                        str(args.GCDType) + '_' + str(i) + '.i3.gz')
    infileList.append(infile)

outfile = dataio.I3File('/home/dvir/workFolder/I3Files/improvedReco/'+ str(args.GCDType) + '/NuGen_improvedReco_' + 
                        str(args.GCDType) + '_improvedRecoTest.i3.gz', 'w')
gcdfile = dataio.I3File(gcdPath)
geometry = gcdfile.pop_frame()["I3Geometry"]

comparisonStringList = [i*2 for i in range(0,15)]
comparisonStringList.extend([i*2+4*30 in range(1,15)])


