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

def getNeighbours(omkey):
    