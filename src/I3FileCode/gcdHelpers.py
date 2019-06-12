#!/usr/bin/env python

'''
a few helper functions that will be commonly used when making
GCD files to test prototype geometries
'''

from icecube import dataio, dataclasses, icetray
from icecube.dataclasses import I3Constants
from icecube.icetray import OMKey, I3Units

cdfile = dataio.I3File(
    '/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/generated/gcd/Calib_and_DetStat_File.i3.gz')
cdframe = cdfile.pop_frame()
start_date = cdframe["I3Calibration"].start_date
end_date = cdframe["I3Calibration"].end_date

# Takes the distance from the surface and returns the z position in
# Icecube coordinates.
# @Param:
# depth:        a float representing vertical distance from the surface
# @Return: 
# the corresponding z value in the IceCube Coordinate system
def convertDepthToZ(depth):
    return I3Constants.SurfaceElev - I3Constants.OriginElev - depth

# Creates a frame with all the needed fields in a C frame
# @Param: 
# geometry:     the I3Geometry object inputted to the frame
# @Return:
# An I3Frame object to be used as the C frame
def generateCFrame(geometry):
    # intialize frame as a C frame
    frame = icetray.I3Frame(icetray.I3Frame.Calibration)

    # add key-value pairs
    frame["I3Geometry"] = geometry
    frame["I3Calibration"] = cdframe["I3Calibration"]
    frame["SPEAbove"] = cdframe["SPEAbove"]
    frame["SPEScalingFactors"] = cdframe["SPEScalingFactors"]

    return frame

# Creates a frame with all the needed fields in a D frame
# @Param: 
# geometry:     the I3Geometry object inputted to the frame
# @Return:
# An I3Frame object to be used as the D frame
def generateDFrame(geometry):
    # intialize frame as a D frame
    frame = icetray.I3Frame(icetray.I3Frame.DetectorStatus)

    # add key-value pairs
    frame["I3Geometry"] = geometry
    frame["I3Calibration"] = cdframe["I3Calibration"]
    frame["SPEAbove"] = cdframe["SPEAbove"]
    frame["SPEScalingFactors"] = cdframe["SPEScalingFactors"]
    frame["I3DetectorStatus"] = cdframe["I3DetectorStatus"]
    frame["BadDomsList"] = cdframe["BadDomsList"]
    frame["BadDomsListSLC"] = cdframe["BadDomsListSLC"]

    return frame

# Generates a string of DOMs in a specified direction. 
# This is represented by an I3OMGeoMap object, which 
# has an OMKey object for a DOM and an I3OMGeo object 
# to represent its geometry.
# @Param: 
# stringNumber: an integer representing the string id
# startPos:     an I3Position object with the location 
#               of the top of the string
# numDoms:      number of DOMs on the string
# spacing:      vertical distance between adjacent DOMs
# direction:    an I3Direction object representing the 
#               direction of the line of the string, 
#               facing away from 
# @Return:
# an I3OMGeoMap object with the DOMs on the string and their
# geometries
def generateOMString( stringNumber, startPos, numDoms, spacing, direction):
    orientation = dataclasses.I3Orientation(0, 0, -1, 1, 0, 0)          # same orientation as icecube DOMs (dir=down)
    area = 0.04439999908208847*I3Units.meter2                           # same area as icecube DOMs
    geomap = dataclasses.I3OMGeoMap()
    x = startPos.x
    y = startPos.y
    z = startPos.z
    dx = spacing*direction.x
    dy = spacing*direction.y
    dz = spacing*direction.z

    # create OMKeys and I3OMGeo for DOMs on string and add them to the map
    for i in xrange(0, numDoms):
        omkey = OMKey(stringNumber, i, 0)
        omGeometry = dataclasses.I3OMGeo()
        omGeometry.omtype = dataclasses.I3OMGeo.OMType.IceCube
        omGeometry.orientation = orientation
        omGeometry.area = area
        omGeometry.position = dataclasses.I3Position(x + dx*i, y + dy*i, z + dz*i)
        geomap[omkey] = omGeometry
    
    return geomap