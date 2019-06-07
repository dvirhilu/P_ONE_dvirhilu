#!/usr/bin/env python

from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import argparse

parser = argparse.ArgumentParser(description = "Generate a simple detector geometry")
parser.add_argument('-l', '--islocal', dest = 'isLocal', 
                    default = 't', help = "configures paths depending on if code is running locally or in cedar (t or f)" )
                    
if args.isLocal == 't':
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz')
else:
	infile = dataio.I3File('/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz')

infile.pop_frame()
infile.pop_frame()

i3Calibration = infile.pop_frame()["I3Calibration"]
dom_cal = i3Calibration.dom_cal
vem_cal = i3Calibration.vem_cal

allI3DOMCalibrations = dom_cal.values()
allI3VEMCalibrations = vem_cal.values()

newDOMCalib = dataclasses.I3DOMCalibration()
# cascadia basin temperature was 
newDOMCalib.temperature = 275 * I3Units.kelvin
newDOMCalib.fadc_gain = sum([domcal.fadc_gain for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.fadc_delta_t = sum([domcal.fadc_delta_t for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.fadc_gain = sum([domcal.beacon_baseline for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.fadc_gain = sum([domcal.fadc_gain for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
