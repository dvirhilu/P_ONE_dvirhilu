#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "Find neutrino distribution from NuGen simulation")
parser.add_argument( '-n', '--minFileNum', help = "smallest file number used" )
parser.add_argument( '-N', '--maxFileNum', help = "largest file number used")
parser.add_argument( '-r', '--runType', help = "The run type (must align with sim run type)")
args = parser.parse_args()

# open file
infileListStep1 = []
for i in range(int(args.minFileNum), int(args.maxFileNum) + 1):
    infile = dataio.I3File('/home/dvir/workFolder/P_ONE_dvirhilu/I3Files/nugen/nugenStep1/NuGen_step1_' + str(args.runType) + '_000' + str(i) + '.i3.gz')
    infileListStep1.append(infile)
    print('file ' + str(i) + 'done' )

qFrameListStep1 = []
for infile in infileListStep1:
    print("qframessss")
    while infile.more():
        qFrameListStep1.append(infile.pop_daq)

print(len(qFrameListStep1))