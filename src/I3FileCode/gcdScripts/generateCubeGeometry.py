#!/usr/bin/env python

from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import argparse
import gcdHelpers

parser = argparse.ArgumentParser(description = "Generate a simple detector geometry")
parser.add_argument('-d', '--domsPerString', dest = 'domsPerString',
                    default = 10, help = "number of doms in the generated string" )
parser.add_argument('-s', '--basicSpacing', dest = 'basicSpacing',
                    default = 15, help = "Spacing between consecutive DOMs on a string. Spacing between strings" )
parser.add_argument('-z', '--depth', dest = 'depth', 
                    default = 1600, help = "starting z position for the geometry" )
args = parser.parse_args()

outfileName = "cubeGeometry_" + str(args.depth) + "_" + str(args.domsPerString) + "_" + str(args.basicSpacing) + ".i3.gz"
outfile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/gcd/cube/' + outfileName, 'w')
domsPerString = int(args.domsPerString)
basicSpacing = float(args.basicSpacing) * I3Units.meter
zpos = gcdHelpers.convertDepthToZ(float(args.depth))
xpos = -0.5*(domsPerString-1)*basicSpacing
ypos = -0.5*(domsPerString-1)*basicSpacing
startingPosition = dataclasses.I3Position(xpos, ypos, zpos*I3Units.meter)

def generateCubeGeometry( domsPerString, startingPos, basicSpacing):
    x = startingPos.x
    y = startingPos.y
    z = startingPos.z
    stringNum = 1
    geomap = dataclasses.I3OMGeoMap()
    stringDirection = dataclasses.I3Direction(0,0,-1)
    spacing = gcdHelpers.generateSpacingList(gcdHelpers.SpacingType.SimpleSpacing, basicSpacing, domsPerString)
    for i in xrange(0, domsPerString):
        for j in xrange( 0, domsPerString):
            startPos = dataclasses.I3Position( x + basicSpacing*i, y + basicSpacing*j, z)
            stringGeoMap = gcdHelpers.generateOMString( stringNum, startPos, domsPerString, spacing, stringDirection)
            geomap.update(stringGeoMap)
            stringNum += 1
    
    return geomap

# generate new geometry object
geometry = dataclasses.I3Geometry()

# fill new geometry
geometry.start_time = gcdHelpers.start_time
geometry.end_time = gcdHelpers.end_time
geometry.omgeo = generateCubeGeometry(domsPerString, startingPosition, basicSpacing)

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