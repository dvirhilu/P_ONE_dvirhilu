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

outfilepath = "/home/dvir/workFolder/P_ONE_dvirhilu/propagationMediumModel/ANTARES/"

# model data
eda = 0.17
scattering = np.array([265, 119])
absorption = np.array([68.6, 23.5])
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
kappa = str(1.0) + '\t# dust absorption wavelength dependence exponent'

# least squares fit to find alpha
coefficientMatrix = np.column_stack( (np.log(wavelengths/400.0), np.ones( len(wavelengths) )) )
leastSquaresMatrix = np.dot(coefficientMatrix.T, coefficientMatrix)
leastSquaresRS = np.dot(coefficientMatrix, scattering.T)
leastSquaresSolution = np.dot( la.inv(leastSquaresMatrix), leastSquaresRS)
print(coefficientMatrix)
print( la.inv(leastSquaresMatrix))
print(leastSquaresSolution)

x = np.linspace(-0.1,0.2,10)
y = [ leastSquaresSolution[0]*i + leastSquaresSolution[1] for i in x]
plt.plot(x,y)
plt.scatter( np.log(wavelengths/400.0), np.log(scattering/400.0) )
plt.show()

# writing data: assumes all absorption falls under pure absorption
min_depth = 1000
max_depth = 2800
spacing = 10
num_elements = ((max_depth - min_depth) / spacing) + 1
depth = np.linspace(1000,2800,num_elements)
