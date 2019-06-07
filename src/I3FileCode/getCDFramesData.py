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
infile.pop_frame()
dframe = infile.pop_frame()

i3Calibration = dframe["I3Calibration"]
i3DetectorStatus = dframe["I3DetectorStatus"]
speabove = dframe["SPEAbove"]
speScalingFactors = dframe["SPEScalingFactors"]

dom_cal = i3Calibration.dom_cal
vem_cal = i3Calibration.vem_cal
dom_status = i3DetectorStatus.dom_status
trigger_status = i3DetectorStatus.trigger_status

allI3DOMCalibrations = dom_cal.values()
allI3VEMCalibrations = vem_cal.values()
allI3DOMStatuses = dom_status.values()

# setting dom calibration (if no known input value average from data)
newDOMCalib = dataclasses.I3DOMCalibration()
nedDOMCalib.tau_parameters.p0 = sum([domcal.tau_parameters.p0 for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
nedDOMCalib.tau_parameters.p1 = sum([domcal.tau_parameters.p1 for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
nedDOMCalib.tau_parameters.p2 = sum([domcal.tau_parameters.p2 for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
nedDOMCalib.tau_parameters.p3 = sum([domcal.tau_parameters.p3 for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
nedDOMCalib.tau_parameters.p4 = sum([domcal.tau_parameters.p4 for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
nedDOMCalib.tau_parameters.p5 = sum([domcal.tau_parameters.p5 for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
nedDOMCalib.tau_parameters.tau_frac = sum([domcal.tau_parameters.tau_frac for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.temperature = 275 * I3Units.kelvin      # cascadia basin temperature was 2 degrees celsius
newDOMCalib.fadc_gain = sum([domcal.fadc_gain for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
nedDOMCalib.fadc_baseline_fit.intercept = sum([domcal.fadc_baseline_fit.intercept for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
nedDOMCalib.fadc_baseline_fit.slope = sum([domcal.fadc_baseline_fit.slope for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.fadc_beacon_baseline = sum([domcal.fadc_beacon_baseline for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.fadc_delta_t = sum([domcal.fadc_delta_t for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.front_end_impedance  = sum([domcal.front_end_impedance for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.atwd_gain[0] = sum([domcal.atwd_gain[0] for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.atwd_gain[1] = sum([domcal.atwd_gain[1] for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.atwd_gain[2] = sum([domcal.atwd_gain[2] for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.transit_time.intercept = sum([domcal.transit_time.intercept for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.transit_time.slope = sum([domcal.transit_time.slope for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.hv_gain_fit.intercept = sum([domcal.hv_gain_fit.intercept for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.hv_gain_fit.slope = sum([domcal.hv_gain_fit.slope for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.dom_cal_version = allI3DOMCalibrations[0].dom_cal_version
ewDOMCalib.atwd_delta_t[0] = sum([domcal.atwd_delta_t[0] for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.atwd_delta_t[1] = sum([domcal.atwd_delta_t[1] for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.spe_disc_calib.intercept = sum([domcal.spe_disc_calib.intercept for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.spe_disc_calib.slope = sum([domcal.spe_disc_calib.intercept.slope for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
ewDOMCalib.mpe_disc_calib.intercept = sum([domcal.mpe_disc_calib.intercept for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.mpe_disc_calib.slope = sum([domcal.mpe_disc_calib.slope for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.pmt_disc_calib.intercept = sum([domcal.pmt_disc_calib.intercept for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.pmt_disc_calib.slope = sum([domcal.pmt_disc_calib.slope for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.relative_dom_eff = 1        # assume same efficiency as icecube doms
newDOMCalib.dom_noise_rate = sum([domcal.dom_noise_rate for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
nedDOMCalib.dom_noise_thermal_rate = sum([domcal.dom_noise_thermal_rate for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
nedDOMCalib.dom_noise_decay_rate = sum([domcal.dom_noise_decay_rate for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.dom_noise_scintillation_mean = sum([domcal.dom_noise_scintillation_mean for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.dom_noise_scintillation_sigma = sum([domcal.dom_noise_scintillation_sigma for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.dom_noise_scintillation_hits  = sum([domcal.dom_noise_scintillation_hits for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.combined_spe_charge_distribution.exp1_amp = sum([domcal.combined_spe_charge_distribution.exp1_amp for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.combined_spe_charge_distribution.exp2_amp = sum([domcal.combined_spe_charge_distribution.exp2_amp for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.combined_spe_charge_distribution.exp1_width = sum([domcal.combined_spe_charge_distribution.exp1_width for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.combined_spe_charge_distribution.exp2_width = sum([domcal.combined_spe_charge_distribution.exp2_width for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.combined_spe_charge_distribution.gaus_amp = sum([domcal.combined_spe_charge_distribution.gaus_amp for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.combined_spe_charge_distribution.gaus_mean = sum([domcal.combined_spe_charge_distribution.gaus_mean for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.combined_spe_charge_distribution.gaus_width = sum([domcal.combined_spe_charge_distribution.gaus_width for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.combined_spe_charge_distribution.compensation_factor = sum([domcal.combined_spe_charge_distribution.compensation_factor for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.combined_spe_charge_distribution.slc_gaus_mean = sum([domcal.combined_spe_charge_distribution.slc_gaus_mean for domcal in allI3DOMCalibrations])/len(allI3DOMCalibrations)
newDOMCalib.combined_spe_charge_distribution.is_valid = True

# setting SPE scaling factor and above
newScalingFactor = sum([scalingfactor for scalingfactor in speScalingFactors.values()])/len(speScalingFactors.values())
newAbove = sum([above for above in speabove.values()])/len(speabove.values())

# setting dom detector status
newDOMStatus = dataclasses.I3DOMStatus()
newDOMStatus.trig_mode = dataclasses.I3DOMStatus.TrigMode.SPE       # standard trigger mode
newDOMStatus.lc_mode = dataclasses.I3DOMStatus.LCMode.UpOrDown      # detect local coincidence from either up or down
newDOMStatus.tx_mode = dataclasses.I3DOMStatus.LCMode.UpAndDown     # send local coincidence info both up and down
newDOMStatus.lc_span = 2                                            # standard local coincidence span
newDOMStatus.lc_window_pre = 1000                                   # value shared across all doms in input file
newDOMStatus.lc_window_post = 1000                                  # value shared across all doms in input file
newDOMStatus.pmt_hv = sum([domstat.pmt_hv for domstat in allI3DOMStatuses])/len(allI3DOMStatuses)
newDOMStatus.spe_threshold = sum([domstat.spe_threshold for domstat in allI3DOMStatuses])/len(allI3DOMStatuses)
newDOMStatus.dac_trigger_bias_0 = sum([domstat.dac_trigger_bias_0 for domstat in allI3DOMStatuses])/len(allI3DOMStatuses)
newDOMStatus.dac_trigger_bias_1 = sum([domstat.dac_trigger_bais_1 for domstat in allI3DOMStatuses])/len(allI3DOMStatuses)
newDOMStatus.dac_trigger_bias_1 = sum([domstat.dac_trigger_bais_1 for domstat in allI3DOMStatuses])/len(allI3DOMStatuses)
newDOMStatus.dac_fadc_ref = 800                                     # value shared across all doms in input file
newDOMStatus.mpe_threshold = sum([domstat.mpe_threshold for domstat in allI3DOMStatuses])/len(allI3DOMStatuses)
newDOMStatus.status_atwd_a = dataclasses.I3DOMStatus.OnOff.On       # allow awtd a readings
newDOMStatus.status_atwd_b = dataclasses.I3DOMStatus.OnOff.On       # allow awtd b readings
newDOMStatus.status_fadc = dataclasses.I3DOMStatus.OnOff.On         # allow fadc readings (off in input file???)
newDOMStatus.delta_compress = dataclasses.I3DOMStatus.OnOff.On      # allow delta compression
newDOMStatus.dom_gain_type = dataclasses.I3DOMStatus.DOMGain.HIGH   # set gain to high
newDOMStatus.slc_active = True                                      # activate SLC readouts
newDOMStatus.fe_pedestal = 2130                                     # value shared across all doms in input file

dummyOMKey = icetray.OMKey(0,0,0)
start_time = dataclasses.I3Time(2019, 0)
end_time = dataclasses.I3Time(2039, 0)
newI3Calibration = dataclasses.I3Calibration()
newI3DetectorStatus = dataclasses.I3DetectorStatus()
newSPEScalingFactors = dataclasses.I3MapKeyDouble()
newSPEAbove = dataclasses.I3MapKeyDouble()
newI3Calibration.dom_cal = dataclasses.Map_OMKey_I3DOMCalibration()
newI3DetectorStatus.dom_status = dataclasses.Map_OMKey_I3DOMStatus()

newI3Calibration.dom_cal[dummyOMKey] = newDOMCalib
