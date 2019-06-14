#!/usr/bin/env python

from optparse import OptionParser
from os.path import expandvars
import os, sys, random



usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="./test_output.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile",default="./test_input.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("-r", "--runnumber", type="string", default="1",
                  dest="RUNNUMBER", help="The run number for this simulation, is used as seed for random generator")
parser.add_option("-f", "--filenr",type="string",default="1",
                  dest="FILENR", help="File number, stream of I3SPRNGRandomService")
parser.add_option("-g", "--gcdfile", default=os.getenv('GCDfile'),
		          dest="GCDFILE", help="Read in GCD file")
parser.add_option("-e","--efficiency", type="float",default=1.2,
                  dest="EFFICIENCY",help="DOM Efficiency ... the same as UnshadowedFraction")
parser.add_option("-n","--noise", default="vuvuzela",
                  dest="NOISE",help="Noise model (vuvuzela/poisson)")
#parser.add_option("-l", "--holeice",  default = "as.flasher_p1_0.30_p2_-1",
#                  dest="HOLEICE", 
#                  help="Pick the hole ice parameterization, corresponds to a file name in $I3_SRC/ice-models/resources/models/angsens/") 
parser.add_option("-m","--icemodel", default="spice_3.2.1",
                  dest="ICEMODEL",help="Ice model (spice_mie, spice_lea, etc)")
# parser.add_option("-s","--scalehad", type="float", default=1.,
#                   dest="SCALEHAD",help="Scale light from hadrons") # This is currently not used

(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

options.FILENR=int(options.FILENR)
options.RUNNUMBER=int(options.RUNNUMBER)

from I3Tray import *
import random

from icecube import icetray, dataclasses, dataio, simclasses #, recclasses
from icecube import phys_services, sim_services, DOMLauncher, DomTools, genie_icetray, clsim, trigger_sim

from RemoveLatePhotons_V5 import RemoveLatePhotons

def BasicHitFilter(frame):
    hits = 0
    if frame.Has("MCPESeriesMap"):
       hits = len(frame.Get("MCPESeriesMap"))
    if hits>0:
#        print "has photons"
        return True
    else:
#       print "does NOT photons"
        return False

def BasicDOMFilter(frame):
    if frame.Has("InIceRawData"):
        if len(frame['InIceRawData']) > 0:
            return True
        else:
            return False
    else:
       return False


tray = I3Tray()

print 'Using RUNNR: ', options.RUNNUMBER
print "DOM efficiency: ", options.EFFICIENCY
print "Using hole ice: ", options.HOLEICE 
print "Looking for ice model in ", expandvars("$I3_SRC/ice-models/resources/models/")

if options.ICEMODEL=="" or options.ICEMODEL==None:
    print "\033[93mNo ice model provided. The baseline efficiency can be found in cfg.txt"
    print "of the ice model used for photon propagation. \033[0m"
    print "\033[93m\033[1mBe very careful! \033[0m"
    icemodel_path=None
else:
    icemodel_path=expandvars("$I3_SRC/ice-models/resources/models/%s"%options.ICEMODEL)
    if os.path.isdir(icemodel_path) :
        print "Folder with ice model found: ", icemodel_path
    else: 
        print "Error! No ice model with such name found in :" 
        print expandvars("$I3_SRC/ice-models/resources/models/")
        exit() 
print "Ice model path: ", icemodel_path
 
# Random service
from globals import max_num_files_per_dataset
tray.AddService("I3SPRNGRandomServiceFactory","sprngrandom")(
    ("Seed",options.RUNNUMBER),
    ("StreamNum",options.FILENR),
    ("NStreams", max_num_files_per_dataset),
    ("instatefile",""),
    ("outstatefile",""),
)

# Now fire up the random number generator with that seed
randomService = phys_services.I3SPRNGRandomService(
    seed = options.RUNNUMBER,
    nstreams = max_num_files_per_dataset,
    streamnum = options.FILENR)


### START ###

tray.AddModule('I3Reader', 'reader',
            FilenameList = [options.GCDFILE, options.INFILE]
            )

 ####
## Remove photons from neutron decay and other processes that take too long (unimportant)
####

tray.AddModule(RemoveLatePhotons, "RemovePhotons",
               InputPhotonSeries = "I3Photons",
               TimeLimit = 1E5) #nanoseconds

####
## Make hits from photons (change efficiency here already!)
####

tray.AddModule("I3GeometryDecomposer", "I3ModuleGeoMap")
# from ReduceHadronicLightyield import HadLightyield

# print "Scaling hadrons with: ", options.SCALEHAD
# tray.AddModule(HadLightyield , "scalecascade", 
# 				Lightyield = options.SCALEHAD)
gcd_file = dataio.I3File(options.GCDFILE)

tray.AddSegment(clsim.I3CLSimMakeHitsFromPhotons, "makeHitsFromPhotons",
#                MCTreeName="I3MCTree_clsim",
#                PhotonSeriesName="UnweightedPhotons2",
                PhotonSeriesName="I3Photons",
                MCPESeriesName="MCPESeriesMap",
                RandomService=randomService,
                DOMOversizeFactor=1.,
                UnshadowedFraction=options.EFFICIENCY,
                IceModelLocation = icemodel_path,
#               UseHoleIceParameterization=holeice
#               HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/%s"%options.HOLEICE),
                GCDFile=gcd_file
                )



#from icecube.BadDomList import bad_dom_list_static
txtfile = os.path.expandvars('$I3_SRC') + '/BadDomList/resources/scripts/bad_data_producing_doms_list.txt'
#BadDoms = bad_dom_list_static.IC86_bad_data_producing_dom_list(118175, txtfile)
tray.AddModule(BasicHitFilter, 'FilterNullMCPE', Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics])
#print BadDoms
mcpe_to_pmt = "MCPESeriesMap"
if options.NOISE == 'poisson':
    print "Error! Poisson noise is not supported anymore. Exiting"
    exit 
elif options.NOISE == 'vuvuzela':
	from icecube import vuvuzela
  # Have removed the vuvuzela traysegment for now, using the module instead below so that the noise parameters are properly taken into account
# 	tray.AddSegment(vuvuzela.AddNoise, 'VuvuzelaNoise',
#                 InputName = mcpe_to_pmt,
#                 OutputName = mcpe_to_pmt + "_withNoise",
# #                ExcludeList = BadDoms,
#                 StartTime = -11*I3Units.microsecond,
#                 EndTime   = 11*I3Units.microsecond,
#                 DisableLowDTCutoff = True
#                 )
	tray.AddModule("Vuvuzela", "vuvuzela_noise" ,
		InputHitSeriesMapName  = mcpe_to_pmt,
		OutputHitSeriesMapName = mcpe_to_pmt + "_withNoise",
	  StartWindow            = -11*I3Units.microsecond,
	  EndWindow              = 11*I3Units.microsecond,
	  IceTop                 = False,
	  InIce                  = True,
	  ScaleFactor            = 1.0,
	  DeepCoreScaleFactor    = 1,
	  DOMsToExclude          = [], # This will be cleaned later by DOM launch cleaner
	  RandomService          = "I3RandomService",
	  SimulateNewDOMs        = True,
	  DisableLowDTCutoff     = True,
	  UseIndividual          = True
	)

        mcpeout = mcpe_to_pmt + '_withNoise'
elif options.NOISE == 'none':
        print '\n*******ERROR: Noiseless simulation!!********\n'
        exit()
        mcpeout = mcpe_to_pmt

else:
	print 'Pick a valid noise model!'
	exit()




tray.AddModule("PMTResponseSimulator","rosencrantz",
    Input=mcpeout,  
    Output=mcpeout + "_weighted",
    MergeHits=True,
    )

tray.AddModule("DOMLauncher", "guildenstern",
    Input= mcpeout + "_weighted",
    Output="InIceRawData_unclean",
    UseTabulatedPT=True,
    )

tray.AddModule("I3DOMLaunchCleaning","launchcleaning")(
       ("InIceInput","InIceRawData_unclean"),
       ("InIceOutput","InIceRawData"),
       ("FirstLaunchCleaning",False),
#       ("CleanedKeys",BadDoms)
       )

# Dropping frames without InIceRawData
tray.AddModule(BasicDOMFilter, 'FilterNullInIce', Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics])
###### triggering 
tray.AddModule('Delete', 'delete_triggerHierarchy',
               Keys = ['I3TriggerHierarchy', 'TimeShift', 'CleanIceTopRawData'])

#time_shift_args = { #'I3MCTreeNames': [],
#                    'I3MCPMTResponseMapNames': [],
#                    'I3MCHitSeriesMapNames' : [] }

tray.AddSegment(trigger_sim.TriggerSim, 'trig', 
                gcd_file = gcd_file,
 #               time_shift_args = time_shift_args,
#added in run_id
                run_id=1)

# Not skipping these keys for now (check what gets dropped in the L2)
skipkeys = ["MCPMTResponseMap",
            "MCTimeIncEventID",
            "I3MCTree_clsim"
            "I3Photons"
            "clsim_stats",
            "InIceRawData_unclean",
            ]

tray.AddModule("I3Writer","writer",
               #SkipKeys=skipkeys, All of these get thrown out by the L2 anyways ... keep them? 
               Filename = options.OUTFILE,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo],
              )

tray.AddModule("TrashCan","adios")

tray.Execute()
tray.Finish()
