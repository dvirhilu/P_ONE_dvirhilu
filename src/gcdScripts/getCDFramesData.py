#!/usr/bin/env python

from icecube import dataio, dataclasses, icetray
from icecube.icetray import OMKey, I3Units
import argparse

parser = argparse.ArgumentParser(description = "Generate a simple detector geometry")
parser.add_argument('-l', '--islocal', dest = 'isLocal', 
                    action = 'store_true', help = "configures paths depending on if code is running locally or in cedar (t or f)" )
args = parser.parse_args()
                    
if args.isLocal:
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/IceCube/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz')
    outfile = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/cal_DS_Files/Calib_and_DetStat_File.i3.gz', 'w')
else:
    infile = dataio.I3File('/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz')
    outfile = dataio.I3File('/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/generated/gcd/Calib_and_DetStat_File.i3.gz', 'w')

# checks if a variable is null or not
def is_nan( variable ):
    return variable != variable

infile.pop_frame()
infile.pop_frame()
infile.pop_frame()
dframe = infile.pop_frame()

i3Calibration = dframe["I3Calibration"]
i3DetectorStatus = dframe["I3DetectorStatus"]
speabove = dframe["SPEAbove"]
speScalingFactors = dframe["SPEScalingFactors"]

dom_cal = i3Calibration.dom_cal
dom_status = i3DetectorStatus.dom_status
trigger_status = i3DetectorStatus.trigger_status

allI3DOMCalibrations = dom_cal.values()
allI3DOMStatuses = dom_status.values()

# setting dom calibration (if no known input value average from data).
# some fields had doms that had nan values. Those values were ignored.
# one of the doms was used instead of creating a new one because some
# attributes were immutable and would have stayed nan if they were 
# intializaed in a new object
newDOMCalib = allI3DOMCalibrations[1000]

tauparam = dataclasses.TauParam()
tauparam.p0 = sum([domcal.tau_parameters.p0 for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.tau_parameters.p0)])/len(allI3DOMCalibrations)
tauparam.p1 = sum([domcal.tau_parameters.p1 for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.tau_parameters.p1)])/len(allI3DOMCalibrations)
tauparam.p2 = sum([domcal.tau_parameters.p2 for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.tau_parameters.p2)])/len(allI3DOMCalibrations)
tauparam.p3 = sum([domcal.tau_parameters.p3 for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.tau_parameters.p3)])/len(allI3DOMCalibrations)
tauparam.p4 = sum([domcal.tau_parameters.p4 for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.tau_parameters.p4)])/len(allI3DOMCalibrations)
tauparam.p5 = sum([domcal.tau_parameters.p5 for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.tau_parameters.p5)])/len(allI3DOMCalibrations)
tauparam.tau_frac = sum([domcal.tau_parameters.tau_frac for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.tau_parameters.tau_frac)])/len(allI3DOMCalibrations)
newDOMCalib.tau_parameters = tauparam

newDOMCalib.temperature = 275 * I3Units.kelvin      # cascadia basin temperature was 2 degrees celsius
newDOMCalib.fadc_gain = sum([domcal.fadc_gain for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.fadc_gain)])/len(allI3DOMCalibrations)

fadcBaselineFit = dataclasses.LinearFit()
fadcBaselineFit.intercept = sum([domcal.fadc_baseline_fit.intercept for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.fadc_baseline_fit.intercept)])/len(allI3DOMCalibrations)
fadcBaselineFit.slope = sum([domcal.fadc_baseline_fit.slope for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.fadc_baseline_fit.slope)])/len(allI3DOMCalibrations)
newDOMCalib.fadc_baseline_fit = fadcBaselineFit

newDOMCalib.fadc_beacon_baseline = sum([domcal.fadc_beacon_baseline for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.fadc_beacon_baseline)])/len(allI3DOMCalibrations)
newDOMCalib.fadc_delta_t = sum([domcal.fadc_delta_t for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.fadc_delta_t)])/len(allI3DOMCalibrations)
newDOMCalib.front_end_impedance  = sum([domcal.front_end_impedance for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.front_end_impedance)])/len(allI3DOMCalibrations)
newDOMCalib.atwd_gain[0] = sum([domcal.atwd_gain[0] for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.atwd_gain[0])])/len(allI3DOMCalibrations)
newDOMCalib.atwd_gain[1] = sum([domcal.atwd_gain[1] for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.atwd_gain[1])])/len(allI3DOMCalibrations)
newDOMCalib.atwd_gain[2] = sum([domcal.atwd_gain[2] for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.atwd_gain[2])])/len(allI3DOMCalibrations)

transitTime = dataclasses.LinearFit()
transitTime.intercept = sum([domcal.transit_time.intercept for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.transit_time.intercept)])/len(allI3DOMCalibrations)
transitTime.slope = sum([domcal.transit_time.slope for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.transit_time.slope)])/len(allI3DOMCalibrations)
newDOMCalib.transit_time = transitTime

hvGainFit = dataclasses.LinearFit()
hvGainFit.intercept = sum([domcal.hv_gain_fit.intercept for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.hv_gain_fit.intercept)])/len(allI3DOMCalibrations)
hvGainFit.slope = sum([domcal.hv_gain_fit.slope for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.hv_gain_fit.slope)])/len(allI3DOMCalibrations)
newDOMCalib.hv_gain_fit = hvGainFit

newDOMCalib.dom_cal_version = allI3DOMCalibrations[0].dom_cal_version
newDOMCalib.atwd_delta_t[0] = sum([domcal.atwd_delta_t[0] for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.atwd_delta_t[0])])/len(allI3DOMCalibrations)
newDOMCalib.atwd_delta_t[1] = sum([domcal.atwd_delta_t[1] for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.atwd_delta_t[1])])/len(allI3DOMCalibrations)

speDiscCalib = dataclasses.LinearFit()
speDiscCalib.intercept = sum([domcal.spe_disc_calib.intercept for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.spe_disc_calib.intercept)])/len(allI3DOMCalibrations)
speDiscCalib.slope = sum([domcal.spe_disc_calib.slope for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.spe_disc_calib.slope)])/len(allI3DOMCalibrations)
newDOMCalib.spe_disc_calib = speDiscCalib

mpeDiscCalib = dataclasses.LinearFit()
mpeDiscCalib.intercept = sum([domcal.mpe_disc_calib.intercept for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.mpe_disc_calib.intercept)])/len(allI3DOMCalibrations)
mpeDiscCalib.slope = sum([domcal.mpe_disc_calib.slope for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.mpe_disc_calib.slope)])/len(allI3DOMCalibrations)
newDOMCalib.mpe_disc_calib = mpeDiscCalib

pmtDiscCalib = dataclasses.LinearFit()
pmtDiscCalib.intercept = sum([domcal.pmt_disc_calib.intercept for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.pmt_disc_calib.intercept)])/len(allI3DOMCalibrations)
pmtDiscCalib.slope = sum([domcal.pmt_disc_calib.slope for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.pmt_disc_calib.slope)])/len(allI3DOMCalibrations)
newDOMCalib.pmt_disc_calib = pmtDiscCalib

newDOMCalib.relative_dom_eff = 1        # assume same efficiency as icecube doms
newDOMCalib.dom_noise_rate = sum([domcal.dom_noise_rate for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.dom_noise_rate)])/len(allI3DOMCalibrations)
newDOMCalib.dom_noise_thermal_rate = sum([domcal.dom_noise_thermal_rate for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.dom_noise_thermal_rate)])/len(allI3DOMCalibrations)
newDOMCalib.dom_noise_decay_rate = sum([domcal.dom_noise_decay_rate for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.dom_noise_decay_rate)])/len(allI3DOMCalibrations)
newDOMCalib.dom_noise_scintillation_mean = sum([domcal.dom_noise_scintillation_mean for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.dom_noise_scintillation_mean)])/len(allI3DOMCalibrations)
newDOMCalib.dom_noise_scintillation_sigma = sum([domcal.dom_noise_scintillation_sigma for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.dom_noise_scintillation_sigma)])/len(allI3DOMCalibrations)
newDOMCalib.dom_noise_scintillation_hits  = sum([domcal.dom_noise_scintillation_hits for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.dom_noise_scintillation_hits)])/len(allI3DOMCalibrations)

combinedSPEChargeDistribution = dataclasses.SPEChargeDistribution()
combinedSPEChargeDistribution.exp1_amp = sum([domcal.combined_spe_charge_distribution.exp1_amp for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.combined_spe_charge_distribution.exp1_amp)])/len(allI3DOMCalibrations)
combinedSPEChargeDistribution.exp1_width = sum([domcal.combined_spe_charge_distribution.exp1_width for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.combined_spe_charge_distribution.exp1_width)])/len(allI3DOMCalibrations)
combinedSPEChargeDistribution.exp2_amp = sum([domcal.combined_spe_charge_distribution.exp2_amp for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.combined_spe_charge_distribution.exp2_amp)])/len(allI3DOMCalibrations)
combinedSPEChargeDistribution.exp2_width = sum([domcal.combined_spe_charge_distribution.exp2_width for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.combined_spe_charge_distribution.exp2_width)])/len(allI3DOMCalibrations)
combinedSPEChargeDistribution.gaus_amp = sum([domcal.combined_spe_charge_distribution.gaus_amp for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.combined_spe_charge_distribution.gaus_amp)])/len(allI3DOMCalibrations)
combinedSPEChargeDistribution.gaus_mean = sum([domcal.combined_spe_charge_distribution.gaus_mean for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.combined_spe_charge_distribution.gaus_mean)])/len(allI3DOMCalibrations)
combinedSPEChargeDistribution.gaus_width = sum([domcal.combined_spe_charge_distribution.gaus_width for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.combined_spe_charge_distribution.gaus_width)])/len(allI3DOMCalibrations)
combinedSPEChargeDistribution.compensation_factor = 1 # no data for spe distribution for these DOMs so use 1 for now
combinedSPEChargeDistribution.slc_gaus_mean = sum([domcal.combined_spe_charge_distribution.slc_gaus_mean for domcal in allI3DOMCalibrations 
                    if not is_nan(domcal.combined_spe_charge_distribution.slc_gaus_mean)])/len(allI3DOMCalibrations)
newDOMCalib.combined_spe_charge_distribution = combinedSPEChargeDistribution

freqQuadFit1 = dataclasses.QuadraticFit()
freqQuadFit2 = dataclasses.QuadraticFit()
dataFreqQuadFit1 = [domcal.atwd_freq_fit[0] for domcal in allI3DOMCalibrations]
dataFreqQuadFit2 = [domcal.atwd_freq_fit[1] for domcal in allI3DOMCalibrations]
freqQuadFit1.quad_fit_a = sum([item.quad_fit_a for item in dataFreqQuadFit1 
                    if not is_nan(item.quad_fit_a)])/len(dataFreqQuadFit1)
freqQuadFit1.quad_fit_b = sum([item.quad_fit_b for item in dataFreqQuadFit1 
                    if not is_nan(item.quad_fit_b)])/len(dataFreqQuadFit1)
freqQuadFit1.quad_fit_c = sum([item.quad_fit_c for item in dataFreqQuadFit1 
                    if not is_nan(item.quad_fit_c)])/len(dataFreqQuadFit1)
freqQuadFit2.quad_fit_a = sum([item.quad_fit_a for item in dataFreqQuadFit2 
                    if not is_nan(item.quad_fit_a)])/len(dataFreqQuadFit2)
freqQuadFit2.quad_fit_b = sum([item.quad_fit_b for item in dataFreqQuadFit2 
                    if not is_nan(item.quad_fit_b)])/len(dataFreqQuadFit2)
freqQuadFit2.quad_fit_c = sum([item.quad_fit_c for item in dataFreqQuadFit2 
                    if not is_nan(item.quad_fit_c)])/len(dataFreqQuadFit2)
newDOMCalib.atwd_freq_fit[0] = freqQuadFit1
newDOMCalib.atwd_freq_fit[1] = freqQuadFit2

# unable to access these fields 
# newDOMCalib.atwd_beacon_baseline = allI3DOMCalibrations[6].atwd_beacon_baseline
# newDOMCalib.atwd_bin_calib_slope = allI3DOMCalibrations[6].atwd_bin_calib_slope

# setting SPE scaling factor and above
newScalingFactor = sum([scalingfactor for scalingfactor in speScalingFactors.values() 
                    if not is_nan(scalingfactor)])/len(speScalingFactors.values())
newAbove = sum([above for above in speabove.values() 
                    if not is_nan(above)])/len(speabove.values())

# setting dom detector status
newDOMStatus = dataclasses.I3DOMStatus()
newDOMStatus.trig_mode = dataclasses.I3DOMStatus.TrigMode.SPE       # standard trigger mode
newDOMStatus.lc_mode = dataclasses.I3DOMStatus.LCMode.UpOrDown      # detect local coincidence from either up or down
newDOMStatus.tx_mode = dataclasses.I3DOMStatus.LCMode.UpAndDown     # send local coincidence info both up and down
newDOMStatus.lc_span = 2                                            # standard local coincidence span
newDOMStatus.lc_window_pre = 1000                                   # value shared across all doms in input file
newDOMStatus.lc_window_post = 1000                                  # value shared across all doms in input file
newDOMStatus.pmt_hv = sum([domstat.pmt_hv for domstat in allI3DOMStatuses
                    if not is_nan(domstat.pmt_hv)])/len(allI3DOMStatuses)
newDOMStatus.spe_threshold = sum([domstat.spe_threshold for domstat in allI3DOMStatuses 
                    if not is_nan(domstat.spe_threshold)])/len(allI3DOMStatuses)
newDOMStatus.dac_trigger_bias_0 = sum([domstat.dac_trigger_bias_0 for domstat in allI3DOMStatuses 
                    if not is_nan(domstat.dac_trigger_bias_0)])/len(allI3DOMStatuses)
newDOMStatus.dac_trigger_bias_1 = sum([domstat.dac_trigger_bias_1 for domstat in allI3DOMStatuses 
                    if not is_nan(domstat.dac_trigger_bias_1)])/len(allI3DOMStatuses)
newDOMStatus.dac_fadc_ref = 800                                     # value shared across all doms in input file
newDOMStatus.mpe_threshold = sum([domstat.mpe_threshold for domstat in allI3DOMStatuses 
                    if not is_nan(domstat.mpe_threshold)])/len(allI3DOMStatuses)
newDOMStatus.status_atwd_a = dataclasses.I3DOMStatus.OnOff.On       # allow awtd a readings
newDOMStatus.status_atwd_b = dataclasses.I3DOMStatus.OnOff.On       # allow awtd b readings
newDOMStatus.status_fadc = dataclasses.I3DOMStatus.OnOff.On         # allow fadc readings (off in input file???)
newDOMStatus.delta_compress = dataclasses.I3DOMStatus.OnOff.On      # allow delta compression
newDOMStatus.dom_gain_type = dataclasses.I3DOMStatus.DOMGain.High   # set gain to high
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
newI3Calibration.start_time = start_time
newI3Calibration.end_time = end_time
newI3Calibration.vem_cal = dataclasses.I3VEMCalibrationMap()

newI3DetectorStatus.daq_configuration_name = "P-ONE_Estimate"
newI3DetectorStatus.dom_status[dummyOMKey] = newDOMStatus
newI3DetectorStatus.start_time = start_time
newI3DetectorStatus.end_time = end_time
newI3DetectorStatus.trigger_status = trigger_status

newSPEAbove[dummyOMKey] = newAbove
newSPEScalingFactors[dummyOMKey] = newScalingFactor

frame = icetray.I3Frame(icetray.I3Frame.DetectorStatus)
frame["I3Calibration"] = newI3Calibration
frame["I3DetectorStatus"] = newI3DetectorStatus
frame["SPEAbove"] = newSPEAbove
frame["SPEScalingFactors"] = newSPEScalingFactors
frame["BadDomsList"] = dataclasses.I3VectorOMKey()
frame["BadDomsListSLC"] = dataclasses.I3VectorOMKey()

outfile.push(frame)
outfile.close()