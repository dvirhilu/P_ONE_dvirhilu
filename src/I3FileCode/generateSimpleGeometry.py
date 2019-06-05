#!/usr/bin/env python

from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import argparse

parser = argparse.ArgumentParser(description = "Generate a simple detector geometry")
parser.add_argument('-l', '--islocal', dest = 'isLocal', 
                    default = 'f', help = "configures paths depending on if code is running locally or in cedar (t or f)" )
args = parser.parse_args()

infile = dataio.I3File('/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz')
outfile = dataio.I3File('/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/generated/gcd/simpleGeometry.i3.gz')
if args.isLocal == 't':
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz')
    outfile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/generated/gcd/simpleGeometry.i3.gz')

def generateOMString( stringNumber, topPos, numDoms, spacing ):
    orientation = dataclasses.I3Orientation(0, 0, -1, 1, 0, 0)          # same orientation as icecube DOMs (dir=down)
    area = 0.04439999908208847*I3Units.meter2                           # same area as icecube DOMs
    geomap = I3OMGeoMap()
    x = topPos.x
    y = topPos.y
    z = topPos.z

    # create OMKeys and I3OMGeo for DOMs on string and add them to the map
    for i in xrange(0, numDoms):
        omkey = OMKey(stringNumber, i, 0)
        omgeometry = dataclasses.I3OMGeo()
        omgeometry.orientation = orientation
        omgeometry.area = area
        omgeometry.position = dataclasses.I3Position(x, y, z - spacing*i)
        geomap[omkey] = omgeometry
    
    return geomap



def generateCubeGeometry( sidelength, topLeftCornerPos):
    

# skip I frame
infile.pop_frame()

# get G,C,D frames
gframe = infile.pop_frame()
cframe = infile.pop_frame()
dframe = infile.pop_frame()

# generate new geometry
geometry = dataclasses.I3Geometry()

# set time
geometry.start_time = dataclasses.I3Time(2019, 0)
geometry.end_time = dataclasses.I3Time(2038,10000000000)



