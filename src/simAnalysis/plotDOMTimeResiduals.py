from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib, csv
from simAnalysis import SimAnalysis
from icecube.phys_services import I3Calculator

parser = argparse.ArgumentParser(description = "Allows to check how clean the time hits of the DOMs are by eye")
parser.add_argument( '-H', '--hitThresh', default = 6, help = "threshold of hits for the DOM to be considered")
parser.add_argument( '-D', '--domThresh', default = 6 , help = "threshold of hit DOMs for the frame to be considered")
args = parser.parse_args()

hitThresh = int(args.hitThresh)
domThresh = int(args.domThresh)
'''
infileList = []
for i in range(300,305):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/HorizGeo/NuGen_step3_HorizGeo_' + str(i) + '.i3.gz')
    infileList.append(infile)

gcd = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/corHorizgeo/CorrHorizGeo_n15_b100.0_a18.0_l3_rise_fall_offset_simple_spacing.i3.gz')
geometry = gcd.pop_frame()["I3Geometry"]
domsUsed = [omkey for omkey in geometry.omgeo.keys() if omkey.string < 16 and omkey.om == 1]

for infile in infileList:
    for frame in infile:
        if SimAnalysis.passFrame(frame, domsUsed, hitThresh, domThresh):
            time = []
            x = []
            timehyp = []
    
            frame = SimAnalysis.writeSigHitsMapToFrame(frame, domsUsed, hitThresh, domThresh)
            mcpeMap = frame["MCPESeriesMap_significant_hits"]
            print frame["NuGPrimary"].dir.azimuth / I3Units.deg
        
            for omkey, mcpeList in mcpeMap:
                position = geometry.omgeo[omkey].position
                x.append(position.x)
                timeList = [mcpe.time for mcpe in mcpeList]
                hitTime = min(timeList)

                primary = frame["NuGPrimary"]
                mctree = frame["I3MCTree"]
                muon = dataclasses.I3MCTree.first_child(mctree, primary)
                muon.shape = dataclasses.I3Particle.InfiniteTrack
                actualTime = I3Calculator.cherenkov_time(muon, position)
                timehyp.append(actualTime)

                time.append( I3Calculator.time_residual(muon, position, hitTime) + actualTime )
        
            plt.figure()
            plt.scatter(time, x, label = 'hit time')
            plt.scatter(timehyp, x, label = 'calculated time')
            plt.xlabel("time (ns)")
            plt.ylabel("distance along string (m)")
            plt.title("Hit detection time on single horizontal string")
            plt.legend()


plt.show()
'''

f = open('/home/dvir/workFolder/P_ONE_dvirhilu/propagationMediumModels/MatthewData/STRAWData_SDOM1_violet_20V_2500Hz.csv')
csv_f = csv.reader(f)
csvData = list(csv_f)

# change range to 0-400ns (for normalization)
time = []
hits = []
for i in range(len(csvData[0])):
    if float(csvData[0][i]) >= 0 and float(csvData[0][i]) <= 400:
        time.append(float(csvData[0][i]))
        hits.append(float(csvData[1][i]))
time = np.array(time)
hits = np.array(hits)

# normalize
dt = time[1] - time[0]
hits = hits / (sum(hits)*dt)

infileList = []
for i in range(300,399):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep2/HorizGeo/NuGen_step2_HorizGeo_' + str(i) + '.i3.gz')
    infileList.append(infile)

gcd = dataio.I3File('/home/dvir/workFolder/I3Files/gcd/corHorizgeo/CorrHorizGeo_n15_b100.0_a18.0_l3_rise_fall_offset_simple_spacing.i3.gz')
geometry = gcd.pop_frame()["I3Geometry"]
timeResid = []
photonCount = 0
for infile in infileList:
    for frame in infile:
        if frame.Stop == I3Frame.DAQ:
            photonMap = frame["I3Photons"]

            for modkey, photonList in photonMap:
                omkey = OMKey(modkey.string, modkey.om, 0)
                for photon in photonList:
                    position = geometry.omgeo[omkey].position
                    hitTime = photon.time
                    primary = frame["NuGPrimary"]
                    mctree = frame["I3MCTree"]
                    muon = dataclasses.I3MCTree.first_child(mctree, primary)
                    muon.shape = dataclasses.I3Particle.InfiniteTrack
                    distance = I3Calculator.cherenkov_distance(muon, position)
                    
                    if photon.wavelength / I3Units.nanometer >= 400 and photon.wavelength / I3Units.nanometer <= 410 and distance > 51 and distance < 55:
                        photonCount += 1
                        timeDiff = I3Calculator.time_residual(muon, position, hitTime)
                        smear = np.random.normal(0, 4)
                        timeResid.append( timeDiff + smear )
                        
bins = np.linspace(0,400,401)
print photonCount
plt.figure()
plt.hist(timeResid, bins = bins,  histtype = 'step', log = True, normed = 1, label = 'Simulation')
plt.step(time, hits, where = 'post', label = 'STRAW data')
plt.xlabel("Time Residual (ns)")
plt.ylabel("Normalized Count")
plt.title("Time Residual Distribution")
plt.legend()
plt.show()
