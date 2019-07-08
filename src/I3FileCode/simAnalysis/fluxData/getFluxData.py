#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray, NuFlux
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "NuFlux module not on local build, so find flux using a script run on cedar")
parser.add_argument('-i', '--infile', help = "input file" )
args = parser.parse_args()

# open file
infile = dataio.I3File(args.infile)

# flux model
flux = NuFlux.makeFlux('honda2006').getFlux

# get all Q frames
qframes = []
while infile.more():
    qframes.append(infile.pop_daq())

# make output file
infilePathStrings = args.infile.split("/")
infileName = infilePathStrings[len(infilePathStrings)-1]
infileAttributes = infileName.split("_")

outname = "fluxData_" + infileAttributes[0] + "_" + infileAttributes[2] + "_" + infileAttributes[3] + ".dat"
outpath = "/project/6008051/dvirhilu/P_ONE_dvirhilu/src/I3FileCode/simAnalysis/fluxData/"
outfile = open(outpath + outname, 'w')

outfile.write("# gives flux multipliers to use when finding the weighting in a simulation for nugen.\n")
outfile.write("# this file is a way to get around the fact that local build does not have the NuFlux module\n")

for frame in qframes:
    eventID = frame["I3EventHeader"].event_id
    primary = frame["NuGPrimary"]
    fluxMult = flux(primary.type, primary.energy, np.cos(primary.dir.zenith))
    outfile.write(str(eventID) + '\t' + str(fluxMult) + '\n')


outfile.close()

