#!/usr/bin/env python

from icecube import dataio, dataclasses, simclasses
from icecube.icetray import I3Units, OMKey
from icecube.clsim import I3CLSimFunctionPolynomial
import numpy as np
import matplotlib.pyplot as plt
from os.path import expandvars


infile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/muongun/muongun_step2/MuonGun_step2_139005_000000.i3.bz2')

# get all Q frames
qframes = []
while(infile.more()):
    qframes.append(infile.pop_daq())

# angularDistribution = GetIceCubeDOMAngularSensitivity()
print(angularDistribution)