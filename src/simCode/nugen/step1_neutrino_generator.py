#!/usr/bin/env python

# script to generate high energy neutrino events for P-ONE. Put together from Neutrino Generator documentation resources

from I3Tray import *
from icecube import icetray, dataclasses, phys_services, sim_services, dataio,  earthmodel_service, neutrino_generator, NuFlux
from icecube.icetray import I3Units
from icecube.dataclasses import I3Particle
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "A scripts to run the neutrino generation simulation step using Neutrino Generator")
# simulation parameters
parser.add_argument('-N', '--runNum', help = "number assigned to this specific run", default = 0 )
parser.add_argument('-n', '--numEvents', help = "number of events produced by the simulation" )
parser.add_argument("-o", "--outfile", help="name and path of output file")
parser.add_argument("-s", "--seed", default=1234567, help="seed for random generator") 
# physics parameters
parser.add_argument('-f', "--flavours", default = "NuMu:NuMuBar", help = "neutrino types to be simulated")
parser.add_argument('-R', "--ratios", default = "1:1", help = "the ratios with which the flavours will be produced" )
parser.add_argument('-E', "--energyLog", default = "2:8", help = "the range of the orders of magnitudes of energies in simulation. 1=GeV")
parser.add_argument('-Z', '--zenithRange', default = "0:180", help = "the range of zenith angles spanned in generation")
parser.add_argument("-p", "--powerLawIndex", default=2.0, help="generation power law index")
# generation surface parameters
parser.add_argument('-x', '--cylinderx', help = "the x coordinate of the center of the injection cylinder")
parser.add_argument('-y', '--cylindery', help = "the y coordinate of the center of the injection cylinder")
parser.add_argument('-z', '--cylinderz', help = "the z coordinate of the center of the injection cylinder")
parser.add_argument('-r', '--cylinderRadius', help = "the radius of the injection cylinder")
parser.add_argument('-l', '--cylinderLength', help = "the length of the injection cylinder")
args = parser.parse_args()

# simulation parameters
runNum = int(args.runNum)
numEvents = int(args.numEvents)
seed = int(args.seed)

# neutrino physics parameters
flavours = args.flavours.split(":")
ratios = [float(ratio) for ratio in args.ratios.split(":")]

eRange = args.energyLog.split(":")
logEMin = float(eRange[0])
logEMax = float(eRange[1])

zenithRange = args.zenithRange.split(":")
zenithMin = float(zenithRange[0]) * I3Units.deg
zenithMax = float(zenithRange[1]) * I3Units.deg

powerLawIndex = float(args.powerLawIndex)

# earth model parameters (has to be set even if simmode is DETECTOR?)
earth = ["PREM_mmc"]
material = ["Standard"]
icecapmodel = "IceSheet"

# injection cylinder parameters
cylinder = [float(args.cylinderRadius), float(args.cylinderLength), float(args.cylinderx), float(args.cylindery), float(args.cylinderz)]

# I3Module to make weighting the distribution easier
class FindEventWeight(icetray.I3Module):

    def __init__(self, context):
        icetray.I3Module.__init__(self,context)
        self.AddParameter("FluxModel", "FluxModel", 'honda2006')
        self.AddParameter("NuTypes", "NeutrinoTypes", ["NuMu", "NuMuBar"])
        self.AddParameter("Ratios", "NeutrinoRatios", [1,1])
        self.AddOutBox("OutBox")
    
    def Configure(self):
        self.fluxModel = self.GetParameter("FluxModel")
        self.types = self.GetParameter("NuTypes")
        self.ratios = self.GetParameter("Ratios")

    def parseTypes(self):
        parsedTypes = []

        for typeString in self.types:
            if typeString == "NuE":
                parsedTypes.append(I3Particle.ParticleType.NuE)
            elif typeString == "NuEBar":
                parsedTypes.append(I3Particle.ParticleType.NuEBar)
            elif typeString == "NuMu":
                parsedTypes.append(I3Particle.ParticleType.NuMu)
            elif typeString == "NuMuBar":
                parsedTypes.append(I3Particle.ParticleType.NuMuBar)
            elif typeString == "NuTau":
                parsedTypes.append(I3Particle.ParticleType.NuTau)
            elif typeString == "NuTauBar":
                parsedTypes.append(I3Particle.ParticleType.NuTauBar)
            else:
                raise RuntimeError("Invalid Neutrino Type: " + typeString)
        
        return parsedTypes


    def getTypeRatioMap(self):
        parsedTypes = self.parseTypes()
        ratioSum = sum(self.ratios)
        normedRatios = [ratio/ratioSum for ratio in self.ratios]

        return {ptype:nRatio for ptype, nRatio in zip(parsedTypes, normedRatios)}

    def DAQ(self, frame):
        # get all necessary data
        primary = frame["NuGPrimary"]
        weightDict = frame["I3MCWeightDict"]
        oneWeight = weightDict["OneWeight"]
        n_events = weightDict["NEvents"]
        ptype = primary.type
        penergy = primary.energy
        ptheta = primary.dir.zenith

        # get flux mult
        flux = NuFlux.makeFlux(self.fluxModel).getFlux
        fluxMult = flux(ptype, penergy, ptheta)

        # get type-ratio dict
        typeRatioMap = self.getTypeRatioMap()

        eventWeight = fluxMult*typeRatioMap[ptype]*oneWeight/n_events
        frame["EventWeight"] = dataclasses.I3Double(eventWeight)
        self.PushFrame(frame)

# initalize tray
tray = I3Tray()

# Random number generator
from globals import max_num_files_per_dataset
randomService = phys_services.I3SPRNGRandomService(
        seed = seed, 
        nstreams = max_num_files_per_dataset,
        streamnum = runNum)


tray.AddModule("I3InfiniteSource","streams",
	       Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
	       EventID=1,
	       IncrementEventID=True)

#
# At least EarthModelService & Steering Service are required
#

tray.AddService("I3EarthModelServiceFactory", "EarthModelService",
                EarthModels = earth,
                MaterialModels = material,
                IceCapType = icecapmodel,
                DetectorDepth = 2600*I3Units.m,
                PathToDataFileDir = "")

tray.AddService("I3NuGSteeringFactory", "steering",
                EarthModelName = "EarthModelService",
                NEvents = numEvents,
                SimMode = "Detector",
                VTXGenMode = "NuGen",
                InjectionMode = "surface",
                CylinderParams = cylinder,
                DoMuonRangeExtension = True,
                UseSimpleScatterForm = True,
                MCTreeName = "I3MCTree_NuGen"
                )

#
# New style configuration
# Primary particle is generated by 
# I3NuGDiffuseSource. By default it
# stores a primary particle with name of 
# NuGPrimary, and if a particle exists with 
# this name in frame, I3NeutrinoGenerator 
# propagates the particle without making 
# a new primary.
# (primary name is configuable)
# You may use I3NuGPointSource either.
#
tray.AddModule("I3NuGDiffuseSource","diffusesource", 
               RandomService = randomService,
               SteeringName = "steering",
               NuTypes = flavours,
               PrimaryTypeRatio = ratios,
               GammaIndex = powerLawIndex,
               EnergyMinLog = logEMin,
               EnergyMaxLog = logEMax,
               ZenithMin = zenithMin,
               ZenithMax = zenithMax,
               AzimuthMin = 0,
               AzimuthMax = 360*I3Units.deg,
               ZenithWeightParam = 1.0,
               AngleSamplingMode = "COS"
              )


tray.AddService("I3NuGInteractionInfoDifferentialFactory", "interaction",
                RandomService = randomService,
                SteeringName = "steering",
                TablesDir = "/project/6008051/dvirhilu/P_ONE_dvirhilu/CrossSectionModels",
                CrossSectionModel = "csms_differential_v1.0"
               )

tray.AddModule("I3NeutrinoGenerator","generator",
                RandomService = randomService,
                SteeringName = "steering",
                InjectorName = "diffusesource",
                InteractionInfoName = "interaction",
#                PropagationWeightMode = "AutoDetect",
                InteractionCCFactor = 1.0,
                InteractionNCFactor = 0.0,
                InteractionGRFactor = 0.0
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
    frame.Put('NEvPerFile', icetray.I3Int(numEvents))

### ADDING PROPAGATOR ###
from icecube import PROPOSAL
from os.path import expandvars

propagators = sim_services.I3ParticleTypePropagatorServiceMap()
	
mediadef=expandvars('$I3_BUILD/PROPOSAL/resources/mediadef')
	
muMinusPropagator = PROPOSAL.I3PropagatorServicePROPOSAL(
		mediadef=mediadef,
		cylinderRadius=1200,
		cylinderHeight=1700,
		type=I3Particle.ParticleType.MuMinus)
muPlusPropagator = PROPOSAL.I3PropagatorServicePROPOSAL(
		mediadef=mediadef,
		cylinderRadius=1200,
		cylinderHeight=1700,
		type=I3Particle.ParticleType.MuPlus)
	
propagators[I3Particle.ParticleType.MuMinus] = muMinusPropagator
propagators[I3Particle.ParticleType.MuPlus] = muPlusPropagator
tray.AddModule('I3PropagatorModule', 'muon_propagator',
		PropagatorServices=propagators,
		RandomService=randomService,
		InputMCTreeName="I3MCTree_NuGen",
        OutputMCTreeName="I3MCTree")

tray.AddModule(FindEventWeight, 'event_weight_finder',
        NuTypes = flavours,
        Ratios = ratios)	

SkipKeys = ["I3MCTree_NuGen"]

tray.AddModule("I3Writer","writer",
 	SkipKeys= SkipKeys,
    Filename = args.outfile)

tray.AddModule("TrashCan", "the can")

tray.Execute()
tray.Finish()
