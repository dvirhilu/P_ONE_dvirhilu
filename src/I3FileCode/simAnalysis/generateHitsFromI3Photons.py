#!/usr/bin/env python

from icecube import dataio, dataclasses, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description = "Takes I3Photons from step2 of the simulations and generates DOM hits")
parser.add_argument('-n', '--runNum', dest = 'runNum', help = "number assigned to this specific run", default = 0 )
parser.add_argument('-g', '--isGenie',dest = 'isGenie', action='store_true', help="is this a simulation done with muongun or genie")
parser.add_argument('-t', '--gcdType',dest = 'GCDType', help = "the type of GCD File used in the simulation")
args = parser.parse_args()

if args.GCDType == 'testString':
    gcdPath = '/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/gcd/testStrings/TestString_n15_b100.0_v50.0_l1_simple_spacing.i3.gz'
elif args.GCDType == 'HorizGeo':
    gcdPath = '/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/gcd/uncorHorizGeo/HorizGeo_n10_b100.0_a90.0_l1_rise_fall_offset_exp_r_spacing.i3.gz'
elif args.GCDType == 'IceCube':
    gcdPath = '/home/dvir/workFolder/I3Files/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz'
else:
    raise RuntimeError("Invalid GCD Type")

if args.isGenie:
    outname = 'genie/customGenHitsGenie/Genie_customGenHits_' + str(args.runNum) + '.i3.gz'
    inPath  = '/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/genie/genie_step2/NuMu/NuMu_C_' + str(args.GCDType) + str(args.runNum) + '.i3.zst'
else:
    outname = 'muongun/customGenHitsMuongun/MuonGun_customGenHits_' + str(args.runNum) + '.i3.gz'
    inPath  = '/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/muongun/muongun_step2/MuonGun_step2_'+ str(args.GCDType) + str(args.runNum) + '.i3.zst'

infile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/muongun/muongun_step2/MuonGun_step2_139005_000000.i3.bz2')
geofile = dataio.I3File(gcdPath)
outfile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/' + outname, 'w')
logfile = open("photonProbabilities.txt",'w')

# get geometry
cframe = geofile.pop_frame(I3Frame.Calibration)
geometry = cframe["I3Geometry"]
geoMap = geometry.omgeo
calibration = cframe["I3Calibration"]
calMap = calibration.dom_cal

# get all Q frames
qframes = []
while(infile.more()):
    qframes.append(infile.pop_daq())

# get files detailing DOM characteristics
inFolder = '/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/'
filenameDomEff = 'icecubeDOMEfficiency.dat'
filenameAngAcc = 'icecubeAngularAcceptance.dat'

domEff = np.loadtxt(inFolder + filenameDomEff, unpack = True)
angAcc = np.loadtxt(inFolder + filenameAngAcc)

def getAngularAcceptanceValue(cos_theta):
    sumVal = 0
    for i in range(len(angAcc)):
        sumVal += cos_theta**i * angAcc[i]
    
    return sumVal

def getDOMAcceptanceValue(wavelength):
    # spacing might not be even just search through
    wlens = domEff[0]
    values = domEff[1]

    if(wavelength < wlens[0]):
        raise RuntimeWarning("wavelength too low, using lowest wavelength in range")
        return values[0]
    elif(wavelength > wlens[len(wlens) - 1]):
        raise RuntimeWarning("wavelength too high, using highest wavelength in range")
        return values[len(wlens) - 1]
    else:
        index = 0
        for i in range(len(wlens)):
            if( wavelength < wlens[i+1]):
                index = i
                break
        fraction = (wavelength - wlens[index])/(wlens[index+1] - wlens[index])

        return values[index] + (values[index+1]-values[index])*fraction

def getSurvivalProbability(photon, omkey):
    if omkey in calMap:
        domcal = calMap[omkey]
        relativeDOMEff = domcal.relative_dom_eff * domcal.combined_spe_charge_distribution.compensation_factor
    else:
        relativeDOMEff = 1

    domGeo = geoMap[omkey]
    domDirection = domGeo.direction
    photonDirection = photon.dir
    dotProduct = photonDirection.x*domDirection.x + photonDirection.y*domDirection.y + photonDirection.z*domDirection.z
    # photon coming in, direction coming out. Sign on dot product should be flipped
    # directions are unit vectors already so cos_theta = -dotProduct (due to sign flip)
    probAngAcc = getAngularAcceptanceValue(-dotProduct)

    probDOMAcc = getDOMAcceptanceValue(photon.wavelength)

    return probAngAcc*probDOMAcc*photon.weight*relativeDOMEff     

def survived(photon,omkey):
    probability = getSurvivalProbability(photon, omkey)
    randomNumber = np.random.uniform()

    if(probability > randomNumber):
        return True

    return False

def generateMCPEList(photons, modkey):
    omkey = OMKey(modkey.string, modkey.om, 0)
    mcpeList = simclasses.I3MCPESeries()
    for photon in photons:
        if survived(photon,omkey):
            mcpe = simclasses.I3MCPE()
            #mcpe.id = dataclasses.I3ParticleID(photon.particleMajorID, photon.particleMinorID)
            mcpe.npe = 1
            mcpe.time = photon.time #TODO: change to corrected time
            mcpeList.append(mcpe)
    
    return mcpeList


# TODO: def getCorrectedTime(photon, omkey):

for frame in qframes:
    photonDOMMap = frame["I3Photons"]
    mcpeMap = simclasses.I3MCPESeriesMap()
    for modkey in photonDOMMap.keys():
        mcpeList = generateMCPEList(photonDOMMap[modkey], modkey)
        omkey = OMKey(modkey.string, modkey.om, 0)
        if len(mcpeList) > 0:
            mcpeMap[omkey] = mcpeList
    
    # only add frame to file if a hit was generated
    if len(mcpeMap.keys()) > 0:
        frame["MCPESeriesMap"] = mcpeMap
        outfile.push(frame)


logfile.close()
outfile.close()

# plot DOM Characteristics

costheta = np.linspace(-1,1,100)
wavelength = np.linspace(260,670,100)

yAng = []
yEff = []
for element in costheta:
    yAng.append(getAngularAcceptanceValue(element))
for element in wavelength:
    yEff.append(getDOMAcceptanceValue(element*I3Units.nanometer))


plt.plot(costheta,yAng)
plt.title("Angular Acceptance of DOM")
plt.xlabel("cos of relative angle")
plt.ylabel("hit probability")
plt.show()

plt.plot(wavelength,yEff)
plt.title("Wavelength Acceptance of DOM")
plt.xlabel("Wavelength (nm)")
plt.ylabel("hit probability")
plt.show()
