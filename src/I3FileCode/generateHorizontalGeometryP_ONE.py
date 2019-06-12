#!/usr/bin/env python

'''
Creates a geometry based on the flat, circular, horizontal
string P-ONE design. 
Parameters that can be modified:
1. number of layers of DOMs (default 1)
2. total angle spanned by detector (Default 90 degrees)
3. DOM spacing on each string (default 100m)
4. Number of DOMS per string (default 10)
5. Total DOMs (default 200)
6. spacing between layers (default 50m)
'''

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units, OMKey
import argparse
import gcdHelpers
import numpy as np

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
parser.add_argument('-h', '--heightSpacing', dest = 'heightSpacing', 
                    default = 50, help = "spacing between layers of DOMs" )         
args = parser.parse_args()

# get parsed arguments
domsPerString = int(args.domsPerString)
spacing = float(args.spacing) * I3Units.meter
zpos = gcdHelpers.convertDepthToZ(float(args.depth)) * I3Units.meter
startPos = dataclasses.I3Position(0,0,zpos)
phi = float(args.angle) * I3Units.deg
layers = int(args.layers)
totalDoms = int(args.totalDoms)
heightSpacing = float(args.heightSpacing)

# calculate parameters from inputs
domsPerLayer = totalDoms/layers
stringsPerLayer = domsPerLayer/domsPerString 
dphi = phi/(stringsPerLayer-1)

if args.outname is not None:
    outname = args.outname
else:
    outname = "HorizGeo_d" + str(domsPerString) +"_s" + str(spacing) + "_a" + str(angle) + "_l" + str(layers) + "_simple"

outfile = dataio.I3File('/home/dvirhilu/workFolder/P_ONE_dvirhilu/I3Files/generated/gcd/' + outname, 'w')

def generateLayer(layerNum):
    stringNumber = layerNum*stringsPerLayer
    x = startPos.x
    y = startPos.y
    z = startpos.z + heightSpacing*layerNum
    layerMap = dataclasses.I3OMGeoMap()
    for i in xrange(0, stringsPerLayer):
        direction = dataclasses.I3Direction(np.cos(i*dphi), np.sin(i*dphi), 0)
        stringStart = dataclasses.I3Position(x + 50*direction.x, y + 50*direction.y, z)
        stringMap = gcdHelpers.generateOMString(stringNumber + i, stringStart, domsPerString, spacing, direction)
        layerMap.update(stringMap)
    return layerMap



    