#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units
import argparse
import gcdHelpers
from gcdHelpers import SpacingType

parser = argparse.ArgumentParser(description = "Generate a vertical string setup to test angular acceptance. Each string represents a line of DOMs in the horizontal geometry")
parser.add_argument('-n', '--domsPerLine', dest = 'domsPerLine',
                    default = 15, help = "number of vertical strings in a generated DOM line" )
parser.add_argument('-b', '--basicSpacing', dest = 'basicSpacing',
                    default = 100, help = "Spacing between consecutive vertical strings on a DOM line" )
parser.add_argument('-l', '--layers', dest = 'layers', 
                    default = 1, help = "number of layers of DOMs in detector i.e. number of DOMs per vertical string" )
parser.add_argument('-v', '--verticalSpacing', dest = 'verticalSpacing', 
                    default = 50, help = "spacing between layers of DOMs" )
parser.add_argument('-s', '--SpacingType', dest = 'spacing_type', type = SpacingType, choices = list(SpacingType), 
                    default = "simple_spacing", help = "type of DOM spacing on DOM lines" )                 
args = parser.parse_args()

# parse arguments and set parameters
domsPerLine = int(args.domsPerLine)
basicSpacing = float(args.basicSpacing) * I3Units.meter
zpos = gcdHelpers.convertDepthToZ(float(2600)) * I3Units.meter
ypos = 0.0
layers = int(args.layers)
verticalSpacing = float(args.verticalSpacing)
spacing_type = args.spacing_type

# create name for output file
outname = "TestString_n" + str(domsPerLine)
outname += "_b" + str(basicSpacing)
outname += "_v" + str(verticalSpacing)
outname += "_l" + str(layers)
outname += "_" + str(spacing_type)
outname += ".i3.gz"

outfile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/gcd/testStrings/' + outname, 'w')

# create new geometry object
geometry = dataclasses.I3Geometry()

# fill new geometry
geometry.start_time = gcdHelpers.start_time
geometry.end_time = gcdHelpers.end_time
geometry.omgeo = dataclasses.I3OMGeoMap()

# create spacing distortion
spacing = gcdHelpers.generateSpacingList(spacing_type, basicSpacing, domsPerLine)

for i in xrange(0, layers):
    stringNum = i + 1
    xpos = verticalSpacing*(i - (layers-1)/2)   # centering position
    startPos = dataclasses.I3Position(xpos, ypos, zpos)
    direction = dataclasses.I3Direction(0, 0, 1)
    stringMap = gcdHelpers.generateOMString(stringNum, startPos, domsPerLine, spacing, direction)
    geometry.omgeo.update(stringMap)

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
    