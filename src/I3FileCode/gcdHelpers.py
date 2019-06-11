#!/usr/bin/env python

from icecube import dataio, dataclasses, icetray
from icecube.dataclasses import I3Constants
from icecube.icetray import OMKey, I3Units

cdfile = dataio.I3File(
    '/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/generated/gcd/Calib_and_DetStat_File.i3.gz')
cdframe = cdfile.pop_frame()

# Takes the distance from the surface and returns the z position in
# Icecube coordinates.
# @Param:
# depth - a float representing vertical distance from the surface
# @Return: 
# the corresponding z value in the IceCube Coordinate system
def convertDepthToZ(depth):
    return I3Constants.SurfaceElev - I3Constants.OriginElev - depth


def generateCFrame(geometry):
    # intialize frame as a C frame
    frame = icetray.I3Frame(icetray.I3Frame.Calibration)

    # add key-value pairs
    frame["I3Geometry"] = geometry
    frame["I3Calibration"] = cdframe["I3Calibration"]
    frame["SPEAbove"] = cdframe["SPEAbove"]
    frame["SPEScalingFactors"] = cdframe["SPEScalingFactors"]

    return frame


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
