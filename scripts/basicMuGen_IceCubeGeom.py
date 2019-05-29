#!/usr/bin/env python

import argparse
from icecube import dataio, dataclasses, icetray, phys_services
from I3Tray import I3Tray
from icecube.icetray import I3Units
from os.path import expandvars
from icecube.MuonGun import load_model, StaticSurfaceInjector, Cylinder, OffsetPowerLaw
from icecube.MuonGun.segments import GenerateBundles

# parse paramaters from command line
parser = argparse.ArgumentParser(description = "Basic Simulation of Neutrino Interactions in IceCube")
parser.add_argument('-o', '--outfile', dest = 'OUTFILE')
parser.add_argument('--gcdfile', dest = 'GCDFILE')
parser.add_argument('--nevents', dest = "NEVENTS")
#parser.add_argument('--minEnergy', dest = "MIN_E")
#parser.add_argument('--maxEnergy', dest = "MAX_E")
#parser.add_argument('--seed', dest = "SEED")
#parser.add_argument('--procnum', dest = "PROCNUM")
#parser.add_argument('--nproc', dest = "NPROC")
args = parser.parse_args()

outfile = dataio.I3File(args.OUTFILE,'w')
gcdfile = dataio.I3File(args.GCDFILE)
nevents = int(args.NEVENTS)
#minE = float(args.MIN_E)
#maxE = float(args.MAX_E)
#seed = int(args.SEED)

tray = I3Tray()
#tray.context['I3RandomService'] = phys_services.I3SPRNGRandomService(1, 10000, 1)

# Use Hoerandel as a template for generating muons
model = load_model('Hoerandel5_atmod12_SIBYLL')
# Generate only single muons, no bundles
model.flux.max_multiplicity = 1
# Center the sampling surface on the barycenter of IC79 strings
surface = Cylinder(1600*I3Units.m, 800*I3Units.m, dataclasses.I3Position(31.25, 19.64, 0))
surface = Cylinder(1600*I3Units.m, 800*I3Units.m)
# Draw energies from an E^-2 power law broken at 1 TeV, from 10 TeV to 10 PeV
spectrum = OffsetPowerLaw(2, 1*I3Units.TeV, 10*I3Units.TeV, 10*I3Units.PeV)
# Set up the generator. This gets stored in a special frame for later reference
generator = StaticSurfaceInjector(surface, model.flux, spectrum, model.radius)

tray.AddSegment(GenerateBundles, 'MuonGenerator', Generator=generator, NEvents=nevents, GCDFile=gcdfile)

tray.Add('Dump')
tray.Execute(10)