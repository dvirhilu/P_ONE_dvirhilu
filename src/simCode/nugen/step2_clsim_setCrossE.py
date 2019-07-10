#!/usr/bin/env python

import argparse
from os.path import expandvars
import os, sys, random

# This script will perform a hybridCLSim propagation.
#
# NOTE: There is no bad_dom_cleaning!!!
#       This you still have to do after the propagation!!!

parser = argparse.ArgumentParser(description = "Makes I3Photons from step1 of the simulations and propagates them to DOMs")
parser.add_argument("-o", "--outfile",default="./test_output.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_argument("-i", "--infile",default="./test_input.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_argument("-r", "--runnumber", default="1",
                  dest="RUNNUMBER", help="The run/dataset number for this simulation, is used as seed for random generator")
parser.add_argument("-l", "--filenr",default="1",
                   dest="FILENR", help="File number, stream of I3SPRNGRandomService")
parser.add_argument("-g", "--gcdfile", default=os.getenv('GCDfile'),
		  dest="GCDFILE", help="Read in GCD file")
parser.add_argument("-e","--efficiency",default=1.2, # Using efficiency > 1 as default so we can support systematics sets
                  dest="EFFICIENCY",help="DOM Efficiency ... the same as UnshadowedFraction")
parser.add_argument("-m","--icemodel", default="spice_3.2.1",
                 dest="ICEMODEL",help="Ice model (spice_mie, spice_lea, etc)")
parser.add_argument("-c","--crossenergy", default=30.0,
                  dest="CROSSENERGY",help="The cross energy where the hybrid clsim approach will be used")
parser.add_argument("-t", action="store_true",  dest="GPU", default=False ,help="Run on GPUs or CPUs")


args = parser.parse_args()

args.FILENR=int(args.FILENR)
args.RUNNUMBER=int(args.RUNNUMBER)
if args.GPU:
        CPU=False
else:
        CPU=True

from I3Tray import *
import random
from icecube import icetray, dataclasses, dataio, simclasses
from icecube import phys_services, sim_services
#from icecube import diplopia
from icecube import clsim

photon_series = "I3Photons"
#print 'CUDA devices: ', args.DEVICE
tray = I3Tray()
print 'Using RUNNUMBER: ', args.RUNNUMBER

# Now fire up the random number generator with that seed
from globals import max_num_files_per_dataset
randomService = phys_services.I3SPRNGRandomService(
    seed = args.RUNNUMBER,
    nstreams = max_num_files_per_dataset,
    streamnum = args.FILENR)

tray.AddService("I3SPRNGRandomServiceFactory","sprngrandom")(
        ("Seed",args.RUNNUMBER),
        ("StreamNum",args.FILENR),
        ("NStreams", max_num_files_per_dataset),
        ("instatefile",""),
        ("outstatefile",""),
)

def BasicHitFilter(frame):
    hits = 0
    if frame.Has(photon_series):
       hits = len(frame.Get(photon_series))
    if hits>0:
       return True
    else:
       return False


### START ###

tray.AddModule('I3Reader', 'reader',
            FilenameList = [args.GCDFILE, args.INFILE]
            )


#tray.AddModule(ModMCTree, "modmctree", mctree="I3MCTree",
#	addhadrons=True
#	)

tray.AddModule("I3GeometryDecomposer", "I3ModuleGeoMap")

icemodel_path = expandvars("$I3_SRC/ice-models/resources/models/" + args.ICEMODEL)
print 'Medium model ', "ANTARES"
print "DOM efficiency: ", args.EFFICIENCY
# Only the photons are made. Still have to convert them to hits!
print "Setting cross energy: " , float(args.CROSSENERGY), "GeV"
#tray.AddSegment(clsim_hybrid.I3CLSimMakePhotons, 'goCLSIM',
print "Using CPUs ", CPU
print "Using GPUs ", args.GPU

gcd_file = dataio.I3File(args.GCDFILE)

tray.AddSegment(clsim.I3CLSimMakePhotons, 'goCLSIM',
                UseCPUs=CPU,
                UseGPUs=args.GPU,
#		UseOnlyDeviceNumber=[0],
#                OpenCLDeviceList=[0],
                MCTreeName="I3MCTree",
                OutputMCTreeName="I3MCTree_clsim",
                FlasherInfoVectName=None,
                MMCTrackListName="MMCTrackList",
                PhotonSeriesName=photon_series,
                ParallelEvents=1000, 
                RandomService=randomService,
                IceModelLocation=icemodel_path,
                #UnWeightedPhotons=True, #turn off optimizations
                UseGeant4=True,
                CrossoverEnergyEM=0.1,
		CrossoverEnergyHadron=float(args.CROSSENERGY),
                StopDetectedPhotons=True,
#                UseHoleIceParameterization=False, # Apply it when making hits!
#                HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.flasher_p1_0.30_p2_-1"),
                DoNotParallelize=False,
                DOMOversizeFactor=1., 
                UnshadowedFraction=args.EFFICIENCY, #normal in IC79 and older CLSim versions was 0.9, now it is 1.0
                GCDFile=gcd_file,
                ExtraArgumentsToI3CLSimModule={
                    #"UseHardcodedDeepCoreSubdetector":True, #may save some GPU memory
                    #"EnableDoubleBuffering":True,
                    "DoublePrecision":False, #will impact performance if true
                    "StatisticsName":"clsim_stats",
                    "IgnoreDOMIDs":[],
                    }
                )


# Tested that all frames go through CLSIM. Removing the ones without any hits to save space.
tray.AddModule(BasicHitFilter, 'FilterNullPhotons', Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

SkipKeys = ["I3MCTree_bak"]

tray.AddModule("I3Writer","writer",
               SkipKeys=SkipKeys,
               Filename = args.OUTFILE,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
              )

tray.AddModule("TrashCan","adios")

tray.Execute()
tray.Finish()

