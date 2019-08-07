#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame
import numpy as np
from simAnalysis import SimAnalysis
import argparse

'''
A fast track reconstruction algorithm that takes into account
the Cherenkov cone structure in its calculation (unlike linefit).
Details of the model can be found here:
https://arxiv.org/pdf/1105.4116.pdf

The reconstruction assumes the particle travels with the speed
of light in a vaccuum: p(t) = q + c*(t-t_0)*u
    q is the position vector of the particle at t=t_0
    u is the direction vector of the particle
    p(t) is the position vector of the particle at time t
    c is the speed of light in a vaccuum

Some more variables definitions for the code:
    z_c is the z component of the point of closest approach (on string)
    d_c is the horizontal distance of the point of closest approach
    t_c is the time the particle reaches its closest approach
    t_gamma is the calculated arrival time of the photon at the DOM
    d_gamma is the calculated distance of the photon to the DOM
    theta_gamma is the calculated angle of inclination of the photon w.r.t. the detector
    L_x is the x position of the DOM
    L_y is the y position of the DOM
    L_z is the z position of the DOM

For convenience, t_0 is chosen to be 0 (linefit implementation uses t_0 = 0)
'''

# constants
c = 2.99792458e8 * I3Units.m / I3Units.second
n = 1.34
t_0 = 0

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

outfile = dataio.I3File('/home/dvir/workFolder/I3Files/linefitReco/'+ str(args.GCDType) + '/NuGen_linefitReco_' + str(args.GCDType) + '_' + str(args.minFileNum) + '_' + str(args.maxFileNum) + '.i3.gz', 'w')
gcdfile = dataio.I3File(gcdPath)
geometry = gcdfile.pop_frame()["I3Geometry"]

def get_z_c(q, u, L_x, L_y):
    # multiplying and I3Position by I3Direction takes their dot product
    numerator = q.z - u.z*(q*u) + u.z*(L_x*u.x + L_y*u.y)
    denominator = 1-u.z**2

    return numerator / denominator

def get_t_c(q, u, L_x, L_y, z_c):
    return t_0 + 1/c * (L_x*u.x + L_y*u.y + z_c*u.z - q*u)

