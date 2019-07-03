#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

nominalModel = np.loadtxt('/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/icecubeAngularAcceptanceNominal.dat')
compModel = np.loadtxt('/home/dvir/workFolder/P_ONE_dvirhilu/DOMCharacteristics/icecubeAngularAcceptance.dat')

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