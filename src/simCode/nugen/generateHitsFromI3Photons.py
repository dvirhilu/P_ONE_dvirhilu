#!/usr/bin/env python

from icecube import dataio, dataclasses, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description = "Takes I3Photons from step2 of the simulations and generates DOM hits")
parser.add_argument('-n', '--runNum',  dest = 'runNum', help = "number assigned to this specific run", default = 0 )
parser.add_argument('-s', '--simType', dest = 'simType', help="which sim tool is used?")
parser.add_argument('-g', '--gcdType', dest = 'GCDType', help = "the type of GCD File used in the simulation")
parser.add_argument('-d', '--domType', dest = 'DOMType', help = "the type of DOM used in the simulation")
parser.add_argument('-H', '--hitThresh', help = "number of total hits required to not cut a frame")
parser.add_argument('-D', '--DOMNumThresh', help = "Number of DOMs with hits required to not cut a frame")
parser.add_argument('-f', '--filePath', help = "path of files will depend on whether code is run locally or not")
args = parser.parse_args()

if args.GCDType == 'testString':
    gcdPath = str(args.filePath) + 'I3Files/gcd/testStrings/HorizTestString_n15_b100.0_v50.0_l1_simple_spacing.i3.gz'
elif args.GCDType == 'HorizGeo':
    gcdPath = str(args.filePath) + 'I3Files/gcd/corHorizgeo/CorrHorizGeo_n15_b100.0_a18.0_l3_rise_fall_offset_simple_spacing.i3.gz'
elif args.GCDType == 'IceCube':
    gcdPath = str(args.filePath) + 'I3Files/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz'
elif args.GCDType == 'cube':
    gcdPath = str(args.filePath) + 'I3Files/gcd/cube/cubeGeometry_1600_15_50.i3.gz'
else:
    raise RuntimeError("Invalid GCD Type")

if args.simType == 'genie':
    outname = 'genie/customGenHitsGenie/Genie_customGenHits_' + str(args.runNum) + '.i3.gz'
    inPath  = str(args.filePath) + 'I3Files/genie/genie_step2/NuMu/NuMu_C_' + str(args.GCDType) + str(args.runNum) + '.i3.zst'
elif args.simType == 'muongun':
    outname = 'muongun/customGenHitsMuongun/MuonGun_customGenHits_' + str(args.runNum) + '.i3.gz'
    inPath  = str(args.filePath) + 'I3Files/muongun/muongun_step2/MuonGun_step2_'+ str(args.GCDType) + str(args.runNum) + '.i3.zst'
elif args.simType == 'nugen':
    outname = 'nugen/nugenStep3/HorizGeo/NuGen_step3_' + str(args.GCDType) + '_' + str(int(args.runNum)+100) + '.i3.gz'
    inPath = str(args.filePath) + 'I3Files/nugen/nugenStep2/HorizGeo/NuGen_step2_' + str(args.GCDType) + '_' + str(args.runNum) + '.i3.gz'
else:
    raise RuntimeError("Invalid Simulation Type")


infile = dataio.I3File(inPath)
geofile = dataio.I3File(gcdPath)
outfile = dataio.I3File(str(args.filePath) + 'I3Files/' + outname, 'w')

# get files detailing DOM characteristics
inFolder = str(args.filePath) + 'P_ONE_dvirhilu/DOMCharacteristics/' + args.DOMType + '/'
filenameDomEff = 'DOMEfficiency.dat'
filenameAngAcc = 'AngularAcceptance.dat'

# get geometry
cframe = geofile.pop_frame(I3Frame.Calibration)
geometry = cframe["I3Geometry"]
geoMap = geometry.omgeo
calibration = cframe["I3Calibration"]
calMap = calibration.dom_cal

def getAngularAcceptanceValue(cos_theta):
    angAcc = np.loadtxt(inFolder + filenameAngAcc)
    
    # check if model exists
    if angAcc.size == 0:
        raise RuntimeError("empty angular Acceptance model")
    # ndarrays of size 1 are not iterable
    elif angAcc.size == 1:
        angAcc = [angAcc]
    
    sumVal = 0
    for i in range(len(angAcc)):
        sumVal += cos_theta**i * angAcc[i]
    
    return sumVal

def getDOMAcceptanceValue(wavelength):
    domEff = np.loadtxt(inFolder + filenameDomEff, unpack = True)
    
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
    photonList = []
    for photon in photons:
        if survived(photon,omkey):
            mcpe = simclasses.I3MCPE()
            #mcpe.id = dataclasses.I3ParticleID(photon.particleMajorID, photon.particleMinorID)
            mcpe.npe = 1
            mcpe.time = photon.time #TODO: change to corrected time
            mcpeList.append(mcpe)
            photonList.append(photon)
    
    return mcpeList, photonList

def passFrame(mcpeMap, domThresh, hitThresh):
    domCount = 0
    for mcpeList in mcpeMap.values():
        if len(mcpeList) >= hitThresh:
            domCount += 1
    
    if domCount < domThresh:
        return False
    
    return True


# TODO: def getCorrectedTime(photon, omkey):

while( infile.more() ):
    frame = infile.pop_daq()
    photonDOMMap = frame["I3Photons"]
    mcpeMap = simclasses.I3MCPESeriesMap()
    succPhotonMap = simclasses.I3CompressedPhotonSeriesMap()
    for modkey in photonDOMMap.keys():
        mcpeList, photonList = generateMCPEList(photonDOMMap[modkey], modkey)
        omkey = OMKey(modkey.string, modkey.om, 0)
        if len(mcpeList) > 0:
            mcpeMap[omkey] = mcpeList
            succPhotonMap[modkey] = photonList
    
    # only add frame to file if a hit was generated
    if passFrame(mcpeMap, int(args.DOMNumThresh), int(args.hitThresh)):
        frame["MCPESeriesMap"] = mcpeMap
        frame["SuccPhotonMap"] = succPhotonMap
        frame.Delete("I3Photons")
        outfile.push(frame)


outfile.close()

'''
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
'''