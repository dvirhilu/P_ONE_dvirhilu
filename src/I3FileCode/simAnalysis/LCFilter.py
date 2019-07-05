#!/usr/bin/env python

from icecube import dataio, dataclasses, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "Takes I3Photons from step2 of the simulations and generates DOM hits")
parser.add_argument('-n', '--runNum', dest = 'runNum', help = "number assigned to this specific run", default = 0 )
parser.add_argument('-g', '--isGenie',dest = 'isGenie', action='store_true', help="is this a simulation done with muongun or genie")
parser.add_argument('-t', '--gcdType',dest = 'GCDType', help = "the type of GCD File used in the simulation")
args = parser.parse_args()

if args.GCDType == 'testString':
    gcdPath = '/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/gcd/testStrings/TestString_n15_b100.0_v50.0_l1_simple_spacing.i3.gz'
elif args.GCDType == 'HorizGeo':
    gcdPath = '/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/gcd/uncorHorizGeo/HorizGeo_n10_b100.0_a90.0_l1_rise_fall_offset_exp_r_spacing.i3.gz'
elif args.GCDType == 'IceCube':
    gcdPath = '/home/dvir/workFolder/I3Files/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz'
else:
    raise RuntimeError("Invalid GCD Type")

if args.isGenie:
    outname = 'genie/customGenHitsGenie/Genie_LCFilter_' + str(args.runNum) + '.i3.gz'
    inPath  = 'genie/customGenHitsGenie/Genie_customGenHits_' + str(args.GCDType) + str(args.runNum) + '.i3.zst'
else:
    outname = 'muongun/customGenHitsMuongun/MuonGun_LCFilter_' + str(args.runNum) + '.i3.gz'
    inPath  = 'muongun/customGenHitsMuongun/Muongun_customGenHits_'+ str(args.GCDType) + str(args.runNum) + '.i3.zst'

infile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/' + inPath)
geofile = dataio.I3File(gcdPath)
outfile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/' + outname, 'w')

frame = infile.pop_frame(I3Frame.Geometry)
geometry = frame["I3Geometry"]


def getNeighbours(omkey, lcRadius):
    neighbours = []
    geoMap = geometry.omgeo
    domgeo = geoMap[omkey]

    for key in geoMap:
        if key != omkey:
            compDOMGeo = geoMap[key]
            pos1 = domgeo.position
            pos2 = compDOMGeo.position
            distancesq = (pos1.x-pos2.x)**2 + (pos1.y-pos2.y)**2 + (pos1.z-pos2.z)**2 

            if distancesq < lcRadius**2:
                neighbours.append(distancesq)

    if len(neighbours) == 0:
        raise RuntimeWarning(str(omkey) + "has no neighbours. LCRadius might be too small")
    
    return neighbours


qframes = []
while infile.more():
    qframes.append(infile.pop_daq())


for frame in qframes:
    mcpeMap = frame["MCPESeriesMap"]
    
