#!/usr/bin/env python

from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import argparse

parser = argparse.ArgumentParser(description = "Generate a simple detector geometry")
parser.add_argument('-l', '--islocal', dest = 'isLocal', 
                    default = 't', help = "configures paths depending on if code is running locally or in cedar (t or f)" )
parser.add_argument('-d', '--domsPerString', dest = 'domsPerString',
                    default = 10, help = "number of doms in the generated string" )
parser.add_argument('-s', '--spacing', dest = 'spacing',
                    default = 15, help = "Spacing between consecutive DOMs on a string. Spacing between strings" )
parser.add_argument('-x', '--xPosition', dest = 'xPosition', 
                    default = 0, help = "starting x position for the geometry" )
parser.add_argument('-y', '--yPosition', dest = 'yPosition', 
                    default = 0, help = "starting y position for the geometry" )
parser.add_argument('-z', '--zPosition', dest = 'zPosition', 
                    default = 0, help = "starting z position for the geometry" )
args = parser.parse_args()

if args.isLocal == 't':
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz')
    outfile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/generated/gcd/simpleGeometry.i3.gz', 'w')
else:
	infile = dataio.I3File('/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz')
	outfile = dataio.I3File('/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/generated/gcd/simpleGeometry.i3.gz', 'w')



domsPerString = int(args.domsPerString)
spacing = float(args.spacing) * I3Units.meter
startingPosition = dataclasses.I3Position( float(args.xPosition)*I3Units.meter, float(args.yPosition)*I3Units.meter, float(args.zPosition)*I3Units.meter)

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

def getDefaultdom_cal(omkeys, inputCalibration):
    domCalibMap = dataclasses.Map_OMKey_I3DOMCalibration()
    domcalib = inputCalibration.dom_cal.values()[0]
    for key in omkeys:
        domCalibMap[key] = domcalib
	return domCalibMap

def getDefaultvem_cal(omkeys, inputCalibration):
    vemCalibMap = dataclasses.Map_OMKey_I3VEMCalibration()
    vemcalib = inputCalibration.vem_cal.values()[0]
    for key in omkeys:
        vemCalibMap[key] = vemcalib
	return vemCalibMap
    
def getDefaultdom_status(omkeys, inputDetectorStatus):
    domStatusMap = dataclasses.Map_OMKey_I3DOMStatus()
    domstatus = inputDetectorStatus.dom_status.values()[0]
    for key in omkeys:
        domStatusMap[key] = domstatus
	return domStatusMap

# skip I frame
infile.pop_frame()

# get G,C,D frames
gframe = infile.pop_frame()
cframe = infile.pop_frame()
dframe = infile.pop_frame()

# generate new geometry, calibration, detectorStatus
geometry = dataclasses.I3Geometry()
calibration = dataclasses.I3Calibration()
detectorStatus = dataclasses.I3DetectorStatus()

# fill new geometry
geometry.start_time = dataclasses.I3Time(2019, 0)
geometry.end_time = dataclasses.I3Time(2038,10000000000)
geometry.omgeo = generateCubeGeometry(domsPerString, startingPosition, spacing)

# set time for new calibration and detectorStatus
calibration.start_time = geometry.start_time
detectorStatus.start_time = geometry.start_time
calibration.end_time = geometry.end_time
detectorStatus.end_time = geometry.end_time

# set default values for calibration and detector status
inputCal = cframe["I3Calibration"]
inputDS = dframe["I3DetectorStatus"]
calibration.dom_cal = getDefaultdom_cal(geometry.omgeo.keys(), inputCal)
calibration.vem_cal = getDefaultvem_cal(geometry.omgeo.keys(), inputCal)
detectorStatus.daq_configuration_name = inputDS.daq_configuration_name
detectorStatus.dom_status = getDefaultdom_status(geometry.omgeo.keys(), inputDS)
detectorStatus.trigger_status = inputDS.trigger_status

# generate new SPEScalingFactor and SPEAbove
scalingFactor = dataclasses.I3MapKeyDouble()
above = dataclasses.I3MapKeyDouble()

for key in geometry.omgeo.keys():
    scalingFactor[key] = 0.75
    above[key] = 0.65

# add keys and values to each frame
gframe.Replace("I3Geometry", geometry)

cframe.Replace("I3Geometry", geometry)
cframe.Replace("I3Calibration", calibration)
cframe.Replace("SPEAbove", above)
cframe.Replace("SPEScalingFactors", scalingFactor)

dframe.Replace("I3Geometry", geometry)
dframe.Replace("I3Calibration", calibration)
dframe.Replace("I3DetectorStatus", detectorStatus)
dframe.Replace("SPEAbove", above)
dframe.Replace("SPEScalingFactors", scalingFactor)
dframe.Replace("BadDomsList", dataclasses.I3VectorOMKey())
dframe.Replace("BadDomsListSLC", dataclasses.I3VectorOMKey())

# push frames onto output file
outfile.push(gframe)
outfile.push(cframe)
outfile.push(dframe)

# close output file
outfile.close()