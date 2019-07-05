#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

nominalModel = np.loadtxt('/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/IceCube/icecubeAngularAcceptanceNominal.dat')
compModel = np.loadtxt('/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/MDOM/AngularAcceptance.dat')

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


plt.plot(costheta, yNom, 'b', label = "nominal model")
plt.plot(costheta, yComp, 'r--', label = "competing model")
plt.title("Comparing Angular Acceptance of Two DOM Models")
plt.xlabel("cos of the relative angle")
plt.ylabel("hit probability")
plt.show()