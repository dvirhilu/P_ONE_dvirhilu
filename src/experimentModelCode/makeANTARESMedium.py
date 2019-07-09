#!/usr/bin/env python

'''
uses properties from "Transmission of light in deep sea water at the
site of the Antares neutrino telescope" Paper (in 
propagationMediumModels/documentation/ANTARTESWaterProperties.pdf)
to create the files clsim takes to construct the propagation medium
'''

import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt

outfilepath = "/home/dvir/workFolder/P_ONE_dvirhilu/propagationMediumModels/ANTARES/"

# model data
eda = 0.17
effScatteringCoeff = np.array([1/265, 1/119])
absorptionCoeff = np.array([1/68.6, 1/23.5])
wavelengths = np.array([473, 375])

# configuration  
header = '# configuration file\n'
domRadius = str(5) + '\t# over-R: DOM radius "oversize" scaling factor\n'               
efficiencyCorr = str(1.0) + '\t# overall DOM efficiency correction\n'
functionShape = str(0.0) + '\t# 0=HG; 1=SAM (scattering function shape parameter)\n'
averageAngle = str((1-eda)*0.924) + '\t# g=<cos(theta)> (formula from ANTARES paper)\n'

configArray = [header, domRadius, efficiencyCorr, functionShape, averageAngle]

outfile = open(outfilepath + "cfg.txt", 'w')

for line in configArray:
    outfile.write(line)

outfile.close()

# parameters
# least squares fit to find alpha and kappa
coefficientMatrix = np.column_stack( (np.log(wavelengths/400.0), np.ones( len(wavelengths) )) )
leastSquaresMatrix = np.matmul(coefficientMatrix.T, coefficientMatrix)
leastSquaresRSScat = np.matmul(coefficientMatrix.T, np.log(effScatteringCoeff.T))
leastSquaresRSAbs = np.matmul(coefficientMatrix.T, np.log(absorptionCoeff.T))
scatLSS = np.matmul( la.inv(leastSquaresMatrix), leastSquaresRSScat)
absLSS = np.matmul( la.inv(leastSquaresMatrix), leastSquaresRSAbs)

# check fit works
x = np.linspace(300,500,201)
y = [ (np.e**scatLSS[1])*( (i/400)**scatLSS[0])  for i in x]
plt.plot(x,y)
plt.scatter( wavelengths, effScatteringCoeff )
plt.show()

# check fit works
plt.figure()
x = np.linspace(300,500,201)
y = [ (np.e**absLSS[1])*( (i/400)**absLSS[0])  for i in x]
plt.plot(x,y)
plt.scatter( wavelengths, absorptionCoeff )
plt.show()

alpha = str(scatLSS[0]) + '\t# scattering wavelength dependence power law exponent\n'
kappa = str(absLSS[0]) + '\t# absorption wavelength dependence power law exponent\n'
b_e_400 = scatLSS[1]
a_e_400 = absLSS[1]

# keep only power law portion of absorption coefficient
A = str(0) + '\t# absorption exponential component multiple - ignored\n'
B = str(0) + '\t# absorption exponential component exponent multiple - ignored\n'

configArray = [alpha, kappa, A, B]

outfile = open(outfilepath + "icemodel.par", 'w')

for line in configArray:
    outfile.write(line)

outfile.close()

# writing data: assumes all absorption falls under power law component
min_depth = 1000
max_depth = 2800
spacing = 10
num_elements = int(((max_depth - min_depth) / spacing)) + 1
depth = np.linspace(1000,2800,num_elements)
b_e = [b_e_400 for i in xrange(num_elements)]
a_e = [a_e_400 for i in xrange(num_elements)]
# assume no temperature difference
dT = [0 for i in xrange(num_elements)]

outfile = open(outfilepath + "icemodel.dat", 'w')

for i in xrange(num_elements):
    line = str(depth[i]) + " " + str(b_e[i]) + " " + str(a_e[i]) + " " + str(dT[i]) + '\n'
    outfile.write(line)

outfile.close()
