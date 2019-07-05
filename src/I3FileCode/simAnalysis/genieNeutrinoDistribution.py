#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "Generates plots for the angular and energy distributions of neutrinos after genie step 1")
parser.add_argument('-n', '--runNum', dest = 'runNum', help = "number assigned to this specific run", default = 900 )
args = parser.parse_args()

# open file
infile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/genie/genie_step1/NuMu_C_' + str(args.runNum) + '_step1.i3.zst')

# get all Q frames
qframes = []
while infile.more():
    qframes.append(infile.pop_daq())
