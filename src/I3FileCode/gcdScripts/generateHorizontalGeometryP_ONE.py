#!/usr/bin/env python

'''
Creates a geometry based on the flat, circular, horizontal
string P-ONE design. 
Parameters that can be modified:
1. number of layers of DOMs (default 1)
2. angular acceptance of string
3. DOM spacing on each string (default 150m)
4. Number of DOMS per string (default 10)
5. spacing between layers (default 50m)
6. different distortion types (detailed in gcdHelpers, default simple_offset, simple_spacing)
'''

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units, OMKey
import argparse
import gcdHelpers
from gcdHelpers import OffsetType, SpacingType
import numpy as np

parser = argparse.ArgumentParser(description = "Generate a geometry based on the flat, circular, horiztal P-ONE design")
parser.add_argument('-n', '--domsPerLine', dest = 'domsPerLine',
                    default = 10, help = "number of vertical strings in a generated DOM line" )
parser.add_argument('-b', '--basicSpacing', dest = 'basicSpacing',
                    default = 150, help = "Spacing between consecutive vertical strings on a DOM line" )
parser.add_argument('-z', '--depth', dest = 'depth', 
                    default = 2600, help = "starting z position for the geometry" )
parser.add_argument('-a', '--angularAcceptance', dest = 'angularAcceptance', 
                    default = 10, help = "angular acceptance of a single DOM line" )
parser.add_argument('-l', '--layers', dest = 'layers', 
                    default = 1, help = "number of layers of DOMs in detector i.e. number of DOMs per vertical string" )
parser.add_argument('-v', '--verticalSpacing', dest = 'verticalSpacing', 
                    default = 50, help = "spacing between layers of DOMs" )
parser.add_argument('-o', '--OffsetType', dest = 'offset_type', type = OffsetType, choices = list(OffsetType), 
                    default = "simple_offset", help = "Offset type to starting DOM line position" )
parser.add_argument('-s', '--SpacingType', dest = 'spacing_type', type = SpacingType, choices = list(SpacingType), 
                    default = "simple_spacing", help = "type of DOM spacing on DOM lines" )                 
args = parser.parse_args()

# get parsed arguments
domsPerLine = int(args.domsPerLine)
basicSpacing = float(args.basicSpacing) * I3Units.meter
zpos = gcdHelpers.convertDepthToZ(float(args.depth)) * I3Units.meter
xpos = 0.0
ypos = 0.0
angularAcceptance = float(args.angularAcceptance) * I3Units.deg
layers = int(args.layers)
verticalSpacing = float(args.verticalSpacing)
offset_type = args.offset_type
spacing_type = args.spacing_type

# calculate parameters from inputs
linesPerLayer = np.ceil( (180-dphi) / (2*dphi) ) + 1    # calculates how many horizontal strings are needed to avoid redundancy
totalDOMs = linesPerLayer * layers                      
dphi = 2 * angularAcceptance                            # to avoid overlap angular spacing needs to be 2*angularAcceptance 


# create name for output file
outname = "CorrHorizGeo_n" + str(domsPerLine)
outname += "_b" + str(basicSpacing)
outname += "_a" + str(dphi/I3Units.deg)
outname += "_l" + str(layers)
outname += "_" + str(offset_type)
outname += "_" + str(spacing_type)
outname += ".i3.gz"

outfile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/generated/gcd/' + outname, 'w')

# create new geometry object
geometry = dataclasses.I3Geometry()

# fill new geometry
geometry.start_time = gcdHelpers.start_time
geometry.end_time = gcdHelpers.end_time
geometry.omgeo = dataclasses.I3OMGeoMap()

# generate distortions
offset = gcdHelpers.generateOffsetList(offset_type, linesPerLayer)
spacing = gcdHelpers.generateSpacingList(spacing_type, basicSpacing, domsPerLine)

# loop to create DOM lines
for i in xrange(0, linesPerLayer):
    direction = dataclasses.I3Direction(np.cos(i*dphi), np.sin(i*dphi), 0)
    lineStart = dataclasses.I3Position(x + offset[i]*direction.x, y + offset[i]*direction.y, z)
    startingStringNum = 1 + i*domsPerLine
    lineMap = generateDOMLine(startStringNum, lineStart, spacing, direction, verticalSpacing, domsPerLine)
    geometry.omgeo.update(lineMap)


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
    