#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units, I3Frame
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib

parser = argparse.ArgumentParser(description = "Find neutrino distribution from NuGen simulation")
parser.add_argument( '-n', '--minFileNum', help = "smallest file number used" )
parser.add_argument( '-N', '--maxFileNum', help = "largest file number used")
parser.add_argument( '-r', '--runType', help = "The run type (must align with GCD type)")
args = parser.parse_args()

# open file
infileList = []
for i in range(int(args.minFileNum), int(args.maxFileNum) + 1):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep1/' + str(args.runType) + '/NuGen_step1_' + str(args.runType) + '_' + str(i) + '.i3.gz')
    infileList.append(infile)

