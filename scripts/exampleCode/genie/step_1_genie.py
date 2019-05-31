#!/usr/bin/env python

# This file has very small modifications from the original example script. It adds MMC for NuMu files.

from optparse import OptionParser
from os.path import expandvars
import os, random

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test_genie.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-l", "--filenr",type="string", default="1",
                   dest="FILENR", help="File number, stream of I3SPRNGRandomService")
parser.add_option("-r", "--runnumber", type="string", default="1",
                  dest="RUNNUMBER", help="The run number for this simulation, also the seed for random number simulations")
parser.add_option("-n", "--numevents", type="string", default="100",
                  dest="NUMEVENTS", help="The number of events per run")
parser.add_option("-f", "--flavor", default="NuMu", 
                  dest="FLAVOR", help="The flavor of the neutrino produced")
parser.add_option("", "--energy-range", default="P", 
                  dest="ENERGYRANGE", help = "A=lowest, B=medium, C=high, P=production, CUSTOM=Must provide E, R, and L")
parser.add_option("", "--elow", type="float", default=9999.9, 
                  dest="ELOW", help = "What is the lower energy bound (only used if energyrange is CUSTOM")
parser.add_option("", "--ehigh", type="float", default=9999.9, 
                  dest="EHIGH", help = "What is the upper energy bound (only used if energyrange is CUSTOM")
parser.add_option("", "--radius", type="float", default=800.0, 
                  dest="RADIUS", help = "What is the lower energy bound (only used if energyrange is CUSTOM")
parser.add_option("", "--length", type="float", default=1600.0, 
                  dest="LENGTH", help = "What is the lower energy bound (only used if energyrange is CUSTOM")
# python step_1_genie.py -o ~/scratch/simulation/L2/test_alex.i3.gz -s 4332 -r 1 -n 1000 -f NuMu

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

options.FILENR=int(options.FILENR)
options.RUNNUMBER=int(options.RUNNUMBER)
options.NUMEVENTS=int(options.NUMEVENTS)

if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

print 'Using RUNNUMBER: ', options.RUNNUMBER
if options.FLAVOR == "NuE":
    en_range = [0.0,0.0]
    spectralIndex = 2.0
    radius = 250.0 ; length = 500.0
    if options.ENERGYRANGE=="A":
        spectralIndex = 2.0
        en_range[0]  = 1.0; en_range[1] = 4.0
    elif options.ENERGYRANGE=="B":
        spectralIndex = 2.0
        en_range[0]  = 4.0; en_range[1] = 12.0
    elif options.ENERGYRANGE=="C":
        spectralIndex = 2.0
        radius = 350.0 ; length = 600.0
        en_range[0]  = 12.0; en_range[1] = 100.0
    elif options.ENERGYRANGE=="D":
        spectralIndex = 2.0
        radius = 550.0 ; length = 1000
        en_range[0]  = 100.0; en_range[1] = 9999.9999
    elif options.ENERGYRANGE=="Z":
        spectralIndex = 1.
        en_range[0]  = 0.5; en_range[1] = 1.0
    elif options.ENERGYRANGE=="CUSTOM":
        en_range[0] = options.ELOW; en_range[1] = options.EHIGH
        radius = options.RADIUS ; length = options.LENGTH
    else:
        print "Unknown option, return!"
        exit()

elif options.FLAVOR == "NuMu":
    en_range = [0.0,0.0]
    spectralIndex = 2.
    if options.ENERGYRANGE=="A":
        en_range[0] = 1.0; en_range[1] = 5.0
        radius = 250.0 ; length = 500.0
    elif options.ENERGYRANGE=="B":
        en_range[0] = 5.0; en_range[1] = 80.0
        radius = 400.0 ; length = 900.0
    elif options.ENERGYRANGE=="C":
        en_range[0] = 80.0; en_range[1] = 1000.0
        radius = 450.0 ; length = 1500.0
    elif options.ENERGYRANGE=="D":
        en_range[0] = 1000.0; en_range[1] = 9999.9999
        radius = 550.0 ; length = 1500.0
    elif options.ENERGYRANGE=="P":
        en_range[0]  = 3.0; en_range[1] = 9999.999999
        radius = 330.0 ; length = 1200.0
        print "Setting P option!!! not divided on ranges"
    elif options.ENERGYRANGE=="Z":
        spectralIndex = 1.
        en_range[0] = 0.5; en_range[1] = 1.0
        radius = 250.0 ; length = 500.0
    elif options.ENERGYRANGE=="CUSTOM":
        en_range[0] = options.ELOW; en_range[1] = options.EHIGH
        radius = options.RADIUS ; length = options.LENGTH
    else:
        print "Unknown option, return!"
        exit()
elif options.FLAVOR == "NuTau":
    en_range = [0.0,0.0]
    spectralIndex = 2.   
    if options.ENERGYRANGE=="A":
        spectralIndex = 2.0
        en_range[0]  = 4.0; en_range[1] = 10.0
        radius = 250.0 ; length = 500.0
    elif options.ENERGYRANGE=="B":
        spectralIndex = 2.00
        en_range[0]  = 10.0; en_range[1] = 50.0
        radius = 350.0 ; length = 600.0
    elif options.ENERGYRANGE=="C":
        spectralIndex = 2.0
        en_range[0]  = 50.0; en_range[1] = 1000.0
        radius = 450.0 ; length = 800.0
    elif options.ENERGYRANGE=="D":
        spectralIndex = 2.0
        en_range[0]  = 1000.0; en_range[1] = 9999.999999  
        radius = 550.0 ; length = 1500.0
    elif options.ENERGYRANGE=="CUSTOM":
        en_range[0] = options.ELOW; en_range[1] = options.EHIGH
        radius = options.RADIUS ; length = options.LENGTH
    else: 
        print "Unknown option! Return!"
        exit()
else:
    en_range = [0.0,0.0]
    spectralIndex = 2.
    en_range[0]  = 3.0; en_range[1] = 9999.999999
print "Radius: ", radius
print "Length: ", length
print "Energy range " + options.ENERGYRANGE+ ": [" + str(en_range[0]) + " , " + str(en_range[1]) + " ] GeV"
    

from I3Tray import *
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services
from icecube import genie_icetray

print "imports complete"

load("libsim-services")


tray = I3Tray()

# Random number generator
from globals import max_num_files_per_dataset
randomService = phys_services.I3SPRNGRandomService(
        seed = options.RUNNUMBER, # +10 for files, which want to fail anyway (usually 1 file per a couple of thousands)
        nstreams = max_num_files_per_dataset,
        streamnum = options.FILENR)


tray.AddModule("I3InfiniteSource","streams",
	       Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
	       Year=2009,
	       DAQTime=158100000000000000,
	       RunNumber=options.FILENR, # Get unique [run,event] num for each event within a dataset
	       EventID=1,
	       IncrementEventID=True)

# Generate the neutrino interactions
tray.AddModule("I3GENIEGenerator","genie_generator",
    RandomService = randomService, 
    GENIEPath= expandvars("$GENIE"),
    SplineFilename = expandvars("/project/6008051/hignight/genie_2_12_8_splines/GENIE_2_12_8_Water_splines.xml"),
    LHAPDFPath = expandvars("$I3_BUILD/genie-icetray/resources/PDFsets"),
    NuEnergyMin = en_range[0]*I3Units.GeV, #3, 195
    NuEnergyMax = en_range[1]*I3Units.GeV,
    PowerLawIndex = spectralIndex, # E^-2.5 spectrum
    GenVolRadius = radius*I3Units.m,
    GenVolLength = length*I3Units.m,
    GenVolDepth = 1950.*I3Units.m, # This option does nothing. Changed it manually in file.
    NeutrinoFlavor = options.FLAVOR, 
    NuFraction = 0.70, # to match lifetime
    MaterialDensity = 0.93*I3Units.g/I3Units.cm3, # ice density
    TargetMixIngredients = [1000080160,1000010010], # O16, H1
    TargetMixQuantities = [1,2], # H2O (O16->1x, H1->2x)
    ForceSingleProbScale = False,
    NEvents = options.NUMEVENTS,
    SystematicNames = ['MaCCRES','MaNCRES','MaCCQE','MaNCEL','MaCOHpi','AhtBY','BhtBY','CV1uBY','CV2uBY'],
    SystematicSteps = [-2,-1,1,2],
    OutputGST=True,
    PositionShift=dataclasses.I3Position(46.29,-34.88,-330.0),  # Changed from dataclasses.I3Position(40,-50,-40), 
    MCTreeName = "I3MCTree_GENIE"
)


# Set up the Driving Time
time = dataclasses.I3Time()
time.set_mod_julian_time(55697, 0, 0)
def DrivingTime( frame ):
	if "DrivingTime" in frame : 
		del frame["DrivingTime"]
	frame.Put("DrivingTime", time )
def NEvents( frame ):
    if "NEvPerFile" in frame:
        del frame['NEvPerFile']
    frame.Put('NEvPerFile', icetray.I3Int(options.NUMEVENTS))

tray.AddModule(DrivingTime, "dt",
	       Streams = [icetray.I3Frame.DAQ] )
tray.AddModule(NEvents, "ne",
	       Streams = [icetray.I3Frame.DAQ] )
if not options.FLAVOR == 'NuMu':
	# Pass the results into an I3MCTree directly
    tray.Add("Rename", Keys = ["I3MCTree_GENIE", "I3MCTree"] )
if options.FLAVOR == 'NuMu':
    ### ADDING PROPAGATOR ###
	from icecube import PROPOSAL, sim_services
	propagators = sim_services.I3ParticleTypePropagatorServiceMap()
	
	mediadef=expandvars('$I3_BUILD/PROPOSAL/resources/mediadef')
	
	muMinusPropagator = PROPOSAL.I3PropagatorServicePROPOSAL(
			mediadef=mediadef,
			cylinderRadius=1200,
			cylinderHeight=1700,
			type=dataclasses.I3Particle.ParticleType.MuMinus)
	muPlusPropagator = PROPOSAL.I3PropagatorServicePROPOSAL(
			mediadef=mediadef,
			cylinderRadius=1200,
			cylinderHeight=1700,
			type=dataclasses.I3Particle.ParticleType.MuPlus)
		
	propagators[dataclasses.I3Particle.ParticleType.MuMinus] = muMinusPropagator
	propagators[dataclasses.I3Particle.ParticleType.MuPlus] = muPlusPropagator
	tray.AddModule('I3PropagatorModule', 'muon_propagator',
			PropagatorServices=propagators,
			RandomService=randomService,
			InputMCTreeName="I3MCTree_GENIE",
            OutputMCTreeName="I3MCTree")	

#######
### Testing: split and pass to root file for checks
#######
# from icecube.tableio import I3TableWriter
# from icecube.rootwriter import I3ROOTTableService
# rootout = options.OUTFILE[:-6] + '.root'
# root_service = I3ROOTTableService(rootout)
# tray.AddModule("I3NullSplitter", "nullsplit")
# tray.AddModule(I3TableWriter,'table_writer',
#                TableService    = [root_service],
#                SubEventStreams = ['nullsplit','in_ice'],
#                BookEverything = True
#                )
# tray.AddModule("Dump", "ShowMeTheMoney")
#
SkipKeys = ["I3MCTree_GENIE"]

tray.AddModule("I3Writer","writer",
 	SkipKeys= SkipKeys,
    Filename = options.OUTFILE)

tray.AddModule("TrashCan", "the can")

tray.Execute()
tray.Finish()
