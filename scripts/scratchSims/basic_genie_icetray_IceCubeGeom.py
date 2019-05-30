#!/usr/bin/env python

import argparse
from icecube import dataio, dataclasses, icetray, phys_services, genie_icetray
from I3Tray import I3Tray
from icecube.icetray import I3Units
import os.path

parser = argparse.ArgumentParser(description = "Basic Simulation of Neutrino Interactions in IceCube")
parser.add_argument('-o', '--outfile', dest = 'OUTFILE')
#parser.add_argument('--gcdfile', dest = 'GCDFILE')
parser.add_argument('--nevents', dest = "NEVENTS")
parser.add_argument('--minEnergy', dest = "MIN_E")
parser.add_argument('--maxEnergy', dest = "MAX_E")
parser.add_argument('--nuFlavour', dest = "NUFLAVOUR")
parser.add_argument('--seed', dest = "SEED")
#parser.add_argument('--procnum', dest = "PROCNUM")
#parser.add_argument('--nproc', dest = "NPROC")
args = parser.parse_args()

outfile = dataio.I3File(args.OUTFILE,'w')
#gcdfile = dataio.I3File(args.GCDFILE)
nevents = int(args.NEVENTS)
minE = float(args.MIN_E)
maxE = float(args.MAX_E)
nuFlavour = str(args.NUFLAVOUR)
seed = int(args.SEED)

tray = I3Tray()
tray.context["I3RandomServices"] = phys_services.I3GSLRandomService(seed)
tray.Add("I3InfiniteSource")
tray.Add("I3GENIEGenerator",
    RandomService = None, # alternatively, this can be None and the I3RandomService can be installed using tray.AddService()
    SplineFilename = os.path.expandvars("$I3_SRC/genie-icetray/resources/splines/splines_water_2.6.4.xml"),
    LHAPDFPath = os.path.expandvars("$I3_SRC/genie-icetray/resources/PDFsets"),
    NuEnergyMin = minE*I3Units.GeV,
    NuEnergyMax = maxE*I3Units.GeV,
    PowerLawIndex = 2., # E^-2 spectrum
    GenVolRadius = 1200.*I3Units.m,
    GenVolLength = 2000.*I3Units.m,
    GenVolDepth = 1950.*I3Units.m,
    NeutrinoFlavor = nuFlavour, # generates neutrinos and anti-neutrinos (1:1)
    MaterialDensity = 0.93*I3Units.g/I3Units.cm3, # ice density
    TargetMixIngredients = [1000080160,1000010010], # O16, H1
    TargetMixQuantities = [1,2], # H2O (O16->1x, H1->2x)
    ForceSingleProbScale = False,
    NEvents = nevents)
tray.Add("I3GENIEResultDictToMCTree")
tray.Add('Dump')
tray.Execute(10) 