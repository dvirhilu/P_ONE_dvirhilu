#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from icecube.icetray import I3Units

nominalModel = np.loadtxt('/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/IceCube/AngularAcceptance.dat')
compModel = np.loadtxt('/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/MDOM/AngularAcceptance.dat')
wavelengthAcc = np.loadtxt('/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/MDOM/DOMEfficiency.dat', unpack = True)

# check if model exists
if nominalModel.size == 0:
    raise RuntimeError("empty nominal model")
# ndarrays of size 1 are not iterable
elif nominalModel.size == 1:
    nominalModel = [nominalModel]

# repeat for comparison model
if compModel.size == 0:
    raise RuntimeError("empty comparison model")
elif compModel.size == 1:
    compModel = [compModel]

costheta = np.linspace(-1,1,100)

yNom = []
yComp = []

for element in costheta:
    sumVal = 0
    for i in range(len(nominalModel)):
        sumVal += element**i * nominalModel[i]
    yNom.append(sumVal)
    sumVal = 0
    for i in range(len(compModel)):
        sumVal += element**i * compModel[i]
    yComp.append(sumVal)


plt.plot(costheta, yNom, 'r--', label = "IceCube Model")
plt.plot(costheta, yComp, 'b', label = "Used Model")
plt.title("DOM Angular Acceptance")
plt.xlabel(r'$cos_{\theta_{rel}}$')
plt.ylabel("Hit Probability")
plt.legend()
plt.show()

plt.plot(wavelengthAcc[0]/I3Units.nanometer, wavelengthAcc[1])
plt.title("DOM Wavelength Acceptance")
plt.xlabel('Wavelength (nm)')
plt.ylabel("Hit Probability")
plt.show()