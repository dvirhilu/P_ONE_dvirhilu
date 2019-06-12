#!/usr/bin/env python

'''
Creates a geometry based on the flat, circular, horizontal
string P-ONE design. 
Parameters that can be modified:
1. number of layers of DOMs
2. total angle spanned by detector
3. DOM spacing on each string
4. Number of DOMS per string
'''

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units, OMKey
import argparse
import gcdHelpers

parser = argparse.ArgumentParser(description = "Generate a geometry based on the flat, circular, horiztal P-ONE design")
parser.add_argument('-d', '--domsPerString', dest = 'domsPerString',
                    default = 10, help = "number of doms in the generated string" )
parser.add_argument('-s', '--spacing', dest = 'spacing',
                    default = 100, help = "Spacing between consecutive DOMs on a string" )
parser.add_argument('-z', '--depth', dest = 'depth', 
                    default = 2600, help = "starting z position for the geometry" )
parser.add_argument('-a', '--angle', dest = 'angle', 
                    default = 90, help = "angle spanned by the detector (degrees)" )
parser.add_argument('-l', '--layers', dest = 'layers', 
                    default = 1, help = "number of layers of DOMs in detector" )
parser.add_argument('-o', '--outname', dest = 'outname', 
                    help = "starting z position for the geometry" )
parser.add_argument('-t', '--totalDoms', dest = 'totalDoms', 
                    default = 200, help = "total number of doms in the detector" )        
args = parser.parse_args()

# get parsed arguments
domsPerString = int(args.domsPerString)
spacing = float(args.spacing) * I3Units.meter
zpos = gcdHelpers.convertDepthToZ(float(args.depth)) * I3Units.meter
xpos = 0
ypos = 0
phi = float(args.angle) * I3Units.deg
layers = int(args.layers)
totalDoms = int(args.totalDoms)

# calculate parameters from inputs
domsPerLayer = totalDoms/layers
stringPerLayer =  

if args.outname is not None:
    outname = args.outname
else:
    outname = "HorizGeo_d" + str(domsPerString) +"_s" + str(spacing) + "_a" + str(angle) + "_l" + str(layers) + "_simple"

outfile = dataio.I3File('/home/dvirhilu/workFolder/P_ONE_dvirhilu/I3Files/generated/gcd/' + outname, 'w')

def generateLayer():

