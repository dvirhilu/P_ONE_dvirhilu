#!/usr/bin/env python

'''
a few helper functions that will be commonly used when making
GCD files to test prototype geometries
'''

from icecube import dataio, dataclasses, icetray
from icecube.dataclasses import I3Constants
from icecube.icetray import OMKey, I3Units
from enum import Enum

cdfile = dataio.I3File(
    '/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/generated/gcd/Calib_and_DetStat_File.i3.gz')
cdframe = cdfile.pop_frame()
calib = cdframe["I3Calibration"]
start_time = calib.start_time
end_time = calib.end_time

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
def generateOMString(stringNumber, startPos, numDoms, spacing, direction):
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


# Generates a list of offsets to the DOM string starting positions. Different
# offset types are described in OffsetType enum class
# @Param:
# offsetType:       an enum representing the different available offset methods
# length:           the length of the resulting list
# @Return:
# a list containing different offset values for the starting positions of the 
# DOM strings
def generateOffsetList(offsetType, length):
    offsetList = []

    if not isinstance(offsetType, OffsetType):
        raise TypeError('direction must be an instance of OffsetType Enum')

    if(offsetType == OffsetType.Simple):
        offsetList = [50 for i in range(0,length)]
    elif(offsetType == OffsetType.ConstantVarying):
        offset = 0
        for i in xrange(0,length):
            if(offset > 100):
                offset = 20
            offsetList.append(offset)
            offset += 20
    
    return offsetList


# an enum class to keep track of different distortion types
# SimpleOffset:         Constant offset of 50m to prevent overlapping (default)
# LinearResetOffset:    Offset starts at and increases by 20m with each string. 
#                       If it reaches 100m, it resets back to 20m.
# LinearRiseFallOffset: Offset starts at 0 and increases by 20m with each string. 
#                       If it reaches 100m, offset starts decreasing by 20m. 
#                       When reaching 20m, starts increasing again.
# OverRSpacing:         Spacing behaves proportional to 1/r, initially starting   
#                       with the base spacing listed in the arguements.
# OverRSquaredSpacing:  Spacing behaves proportional to 1/r^2, initially starting 
#                       with the base spacing listed in the arguements.
class DistortionType(Enum):
    SimpleOffset = "simple_offset"
    LinearResetOffset = "linear_reset_offset"
    LinearRiseFallOffset = "rise_fall_offset"
    OverRSpacing = "1/r_spacing"
    OverRSquaredSpacing = "1/r^2_spacing"

    def __str__(self):
        return self.value