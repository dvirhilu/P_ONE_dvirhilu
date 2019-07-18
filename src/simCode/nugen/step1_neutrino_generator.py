#!/usr/bin/env python

# script to generate high energy neutrino events for P-ONE. Put together from Neutrino Generator documentation resources

from I3Tray import *
from icecube import icetray, dataclasses, phys_services, sim_services, dataio,  earthmodel_service, neutrino_generator, NuFlux
from icecube.icetray import I3Units, I3Frame
from icecube.dataclasses import I3Particle
from icecube.simclasses import I3MMCTrack
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "A scripts to run the neutrino generation simulation step using Neutrino Generator")
# simulation parameters
parser.add_argument('-N', '--runNum', help = "number assigned to this specific run", default = 0 )
parser.add_argument('-n', '--numEvents', help = "number of events produced by the simulation" )
parser.add_argument('-o', '--outfile', help="name and path of output file")
parser.add_argument('-s', '--seed', default=1234567, help="seed for random generator") 
parser.add_argument('-g', '--gcdFile', help = "gcd file used for simulation set")
# physics parameters
parser.add_argument('-f', '--flavours', default = "NuMu:NuMuBar", help = "neutrino types to be simulated")
parser.add_argument('-R', '--ratios', default = "1:1", help = "the ratios with which the flavours will be produced" )
parser.add_argument('-E', '--energyLog', default = "2:8", help = "the range of the orders of magnitudes of energies in simulation. 1=GeV")
parser.add_argument('-a', '--coneAngle', default = 45, help = "the angle spanned by generation")
parser.add_argument('-p', '--powerLawIndex', default=2.0, help="generation power law index")
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

zenithMin = (90 - float(args.coneAngle)) * I3Units.deg
zenithMax = (90 + float(args.coneAngle)) * I3Units.deg
azimuthMin = -float(args.coneAngle) * I3Units.deg
azimuthMax = float(args.coneAngle) * I3Units.deg

powerLawIndex = float(args.powerLawIndex)

# earth model parameters (has to be set even if simmode is DETECTOR?)
earth = ["PREM_mmc"]
material = ["Standard"]
icecapmodel = "IceSheet"

# injection cylinder parameters
cylinder = [float(args.cylinderRadius), float(args.cylinderLength), float(args.cylinderx), float(args.cylindery), float(args.cylinderz)]

# I3Module to make weighting the distribution easier
class ClosestApproachFilter(icetray.I3Module):

    def __init__(self, context):
        icetray.I3Module.__init__(self,context)
        self.AddParameter("GCDFile", "GCDFile")
        self.AddOutBox("OutBox")
    
    def Configure(self):
        gcdFile = self.GetParameter("GCDFile")

        gcd = dataio.I3File(gcdFile)
        geometry = gcd.pop_frame(I3Frame.Geometry)["I3Geometry"]
        self.geoMap = geometry.omgeo

    def getClosestApproachDistance(self, muon):
        muonPos = muon.pos
        muonDir = muon.dir
        mx = muonPos.x
        my = muonPos.y
        mz = muonPos.z
        mex = muonDir.x
        mey = muonDir.y
        mez = muonDir.z

        for domgeo in self.geoMap.values():
            domPos = domgeo.position
            domx = domPos.x
            domy = domPos.y
            domz = domPos.z

            relPosx = domx - mx
            relPosy = domy - my
            relPosz = domz - mz

            projDist = relPosx*mex + relPosy*mey + relPosz*mez

            projx = mx + mex*projDist
            projy = my + mey*projDist
            projz = mz + mez*projDist

            rSquared = (projx - domx)**2 + (projy - domy)**2 + (projz - domz)**2

            return np.sqrt(rSquared)


    def DAQ(self, frame):
        # get all necessary data
        trackList = frame["MMCTrackList"]
        muon = trackList[0].GetI3Particle()
        closestAppDistance = self.getClosestApproachDistance(muon)

        frame["ClosestAppoachDistance"] = dataclasses.I3Double(closestAppDistance)
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
               AzimuthMin = azimuthMin,
               AzimuthMax = azimuthMax,
               ZenithWeightParam = 1.0,
               AngleSamplingMode = "COS"
              )


tray.AddService("I3NuGInteractionInfoDifferentialFactory", "interaction",
                RandomService = randomService,
                SteeringName = "steering",
                TablesDir = "/home/users/dhilu/P_ONE_dvirhilu/CrossSectionModels",
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

tray.AddModule(ClosestApproachFilter, 'closesr_approach_filter',
        GCDFile = args.gcdFile)	

SkipKeys = ["I3MCTree_NuGen"]

tray.AddModule("I3Writer","writer",
 	SkipKeys= SkipKeys,
    Filename = args.outfile)

tray.AddModule("TrashCan", "the can")

tray.Execute()
tray.Finish()
