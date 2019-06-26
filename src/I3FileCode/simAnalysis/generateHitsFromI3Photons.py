#!/usr/bin/env python

from icecube import dataio, dataclasses, simclasses
from icecube.icetray import I3Units, OMKey
import numpy as np
import matplotlib.pyplot as plt
from os.path import expandvars


infile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/muongun/muongun_step2/MuonGun_step2_139005_000000.i3.bz2')

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

# plot distributions
cos_theta = np.linspace(-1,1,100)
probAcc = np.array([])
for i in range(len(cos_theta)):
    sumVal = 0
    for j in range(len(angAcc)):
        sumVal += (cos_theta[i])**j * angAcc[j]
    probAcc = np.append(probAcc, sumVal)

wavelengths = domEff[0] / I3Units.nanometer
efficiency = domEff[1]

plt.plot(cos_theta, probAcc)
plt.title("Angular Acceptance vs. Incidence Angle Theta")
plt.xlabel("cos(theta)")
plt.ylabel("Angular Acceptance")
plt.show()

plt.figure()
plt.plot(wavelengths, efficiency)
plt.title("DOM Efficiency vs. Wavelength")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Efficiency")
plt.show()