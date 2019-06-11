#!/usr/bin/env python

from icecube import dataio, dataclasses, icetray
from icecube.dataclasses import I3Constants
from icecube.icetray import OMKey, I3Units
import argparse

parser = argparse.ArgumentParser(description = "Generate a simple detector geometry")
parser.add_argument('-l', '--islocal', dest = 'isLocal', 
                    action='store_true', help = "configures paths depending on if code is running locally or in cedar (t or f)" )
parser.add_argument('-d', '--domsPerString', dest = 'domsPerString',
                    default = 10, help = "number of doms in the generated string" )
parser.add_argument('-s', '--spacing', dest = 'spacing',
                    default = 15, help = "Spacing between consecutive DOMs on a string. Spacing between strings" )
parser.add_argument('-x', '--xPosition', dest = 'xPosition', 
                    default = 0, help = "starting x position for the geometry" )
parser.add_argument('-y', '--yPosition', dest = 'yPosition', 
                    default = 0, help = "starting y position for the geometry" )
parser.add_argument('-z', '--depth', dest = 'depth', 
                    default = 1000, help = "starting z position for the geometry" )
args = parser.parse_args()

outfileName = "cubeGeometry_" + str(args.depth) + "_" + str(args.domsPerString) + "_" + str(args.spacing) + ".i3.gz"

if args.isLocal:
    outfile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/generated/gcd/' + outfileName, 'w')
    cdfile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/generated/gcd/Calib_and_DetStat_File.i3.gz')
else:
	outfile = dataio.I3File('/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/generated/gcd/' + outfileName, 'w')
        cdfile = dataio.I3File('/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/generated/gcd/Calib_and_DetStat_File.i3.gz')

domsPerString = int(args.domsPerString)
spacing = float(args.spacing) * I3Units.meter
zpos = I3Constants.SurfaceElev - I3Constants.OriginElev - float(args.depth)
startingPosition = dataclasses.I3Position( float(args.xPosition)*I3Units.meter, float(args.yPosition)*I3Units.meter, zpos*I3Units.meter)

def generateOMString( stringNumber, topPos, numDoms, spacing ):
    orientation = dataclasses.I3Orientation(0, 0, -1, 1, 0, 0)          # same orientation as icecube DOMs (dir=down)
    area = 0.04439999908208847*I3Units.meter2                           # same area as icecube DOMs
    geomap = dataclasses.I3OMGeoMap()
    x = topPos.x
    y = topPos.y
    z = topPos.z

    # create OMKeys and I3OMGeo for DOMs on string and add them to the map
    for i in xrange(0, numDoms):
        omkey = OMKey(stringNumber, i, 0)
        omGeometry = dataclasses.I3OMGeo()
        omGeometry.omtype = dataclasses.I3OMGeo.OMType.IceCube
        omGeometry.orientation = orientation
        omGeometry.area = area
        omGeometry.position = dataclasses.I3Position(x, y, z - spacing*i)
        geomap[omkey] = omGeometry
    
    return geomap

def generateCubeGeometry( domsPerString, startingPos, spacing):
    x = startingPos.x * I3Units.meter
    y = startingPos.y * I3Units.meter
    z = startingPos.z * I3Units.meter
    stringNum = 1
    geomap = dataclasses.I3OMGeoMap()
    for i in xrange(0, domsPerString):
        for j in xrange( 0, domsPerString):
            topPos = dataclasses.I3Position( x + spacing*i, y + spacing*j, z)
            stringGeoMap = generateOMString( stringNum, topPos, domsPerString, spacing)
            geomap.update(stringGeoMap)
            stringNum += 1
    
    return geomap

# get C,D frames
cdframe = cdfile.pop_frame()

# generate new GCD Objects
geometry = dataclasses.I3Geometry()
calibration = cdframe["I3Calibration"]
detectorStatus = cdframe["I3DetectorStatus"]
scalingFactor = cdframe["SPEScalingFactors"]
above = cdframe["SPEAbove"]
badDoms = cdframe["BadDomsList"]
badDomsSLC = cdframe["BadDomsListSLC"]

# fill new geometry
geometry.start_time = calibration.start_time
geometry.end_time = calibration.end_time
geometry.omgeo = generateCubeGeometry(domsPerString, startingPosition, spacing)

# generate new frames
gframe = icetray.I3Frame(icetray.I3Frame.Geometry) 
cframe = icetray.I3Frame(icetray.I3Frame.Calibration)
dframe = icetray.I3Frame(icetray.I3Frame.DetectorStatus)

# add keys and values to each frame
gframe["I3Geometry"] = geometry

cframe["I3Geometry"] = geometry
cframe["I3Calibration"] = calibration
cframe["SPEAbove"] = above
cframe["SPEScalingFactors"] = scalingFactor

dframe["I3Geometry"] = geometry
dframe["I3Calibration"] = calibration
dframe["I3DetectorStatus"] = detectorStatus
dframe["SPEAbove"] = above
dframe["SPEScalingFactors"] = scalingFactor
dframe["BadDomsList"] = badDoms
dframe["BadDomsListSLC"] = badDomsSLC

# push frames onto output file
outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

# close output file
outfile.close()