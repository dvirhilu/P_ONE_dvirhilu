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
7.
'''

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units, OMKey
import argparse
import gcdHelpers
from gcdHelpers import DistortionType
import numpy as np

parser = argparse.ArgumentParser(description = "Generate a geometry based on the flat, circular, horiztal P-ONE design")
parser.add_argument('-n', '--domsPerString', dest = 'domsPerString',
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
parser.add_argument('-v', '--verticalSpacing', dest = 'verticalSpacing', 
                    default = 50, help = "spacing between layers of DOMs" )
parser.add_argument('d', '--DistortionType', dest = 'distType', type = DistortionType, choices = list(DistortionType), 
                    default = ["simple_offset"], action='append', help = "allows for a list of distortions to be provided" )         
args = parser.parse_args()

# get parsed arguments
domsPerString = int(args.domsPerString)
spacing = float(args.spacing) * I3Units.meter
zpos = gcdHelpers.convertDepthToZ(float(args.depth)) * I3Units.meter
startPos = dataclasses.I3Position(0,0,zpos)
phi = float(args.angle) * I3Units.deg
layers = int(args.layers)
totalDoms = int(args.totalDoms)
verticalSpacing = float(args.verticalSpacing)
offset_type = args.offset_type

# calculate parameters from inputs
domsPerLayer = totalDoms/layers
stringsPerLayer = domsPerLayer/domsPerString 
dphi = phi/(stringsPerLayer-1)

if args.outname is not None:
    outname = args.outname
else:
    outname = "HorizGeo_d" + str(domsPerString) +"_s" + str(spacing) + "_a" + str(phi/I3Units.deg) + "_l" + str(layers) + "_" + str(offset_type) + ".i3.gz"

outfile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/generated/gcd/' + outname, 'w')

def generateLayer(layerNum):
    stringNumber = layerNum*stringsPerLayer
    x = startPos.x
    y = startPos.y
    z = startPos.z + verticalSpacing*layerNum
    # offset so that first DOMs in each string don't overlap
    offset = gcdHelpers.generateOffsetList(offset_type, stringsPerLayer)
    layerMap = dataclasses.I3OMGeoMap()
    
    for i in xrange(0, stringsPerLayer):
        # tilt every new string by an angle of dphi
        direction = dataclasses.I3Direction(np.cos(i*dphi), np.sin(i*dphi), 0)
        stringStart = dataclasses.I3Position(x + offset[i]*direction.x, y + offset[i]*direction.y, z)
        stringMap = gcdHelpers.generateOMString(stringNumber + i, stringStart, domsPerString, spacing, direction)
        layerMap.update(stringMap)
    
    return layerMap

# create new geometry object
geometry = dataclasses.I3Geometry()

# fill new geometry
geometry.start_time = gcdHelpers.start_time
geometry.end_time = gcdHelpers.end_time
geometry.omgeo = dataclasses.I3OMGeoMap()

for i in xrange(0,layers):
    layerMap = generateLayer(i)
    geometry.omgeo.update(layerMap)

# generate new frames
gframe = icetray.I3Frame(icetray.I3Frame.Geometry) 
cframe = gcdHelpers.generateCFrame(geometry)
dframe = gcdHelpers.generateDFrame(geometry)

# add keys and values to G frame
gframe["I3Geometry"] = geometry

# push frames onto output file
outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

# close output file
outfile.close()
    