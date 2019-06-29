#!/usr/bin/env python

#########################
#  INTRO
#########################
print("Starting python script")
from os.path import expandvars
from argparse import ArgumentParser

import numpy as np

from I3Tray import I3Tray

import icecube
from icecube import icetray, dataio, dataclasses, phys_services, MuonGun
from icecube.icetray import I3Units
from icecube.MuonGun.segments import GenerateBundles
from icecube.sim_services.propagation import get_propagators
icetray.set_log_level(icetray.logging.I3LogLevel.LOG_WARN)

#########################
#  PARSE OPTIONS
#########################

parser = ArgumentParser()
parser.add_argument("-g", "--gcdfile", type=str, help="GCD file path (.i3{.gz} format)")
parser.add_argument("-o", "--outfile", type=str, help="Output file path (.i3{.gz} format)")
#parser.add_argument("-d", "--dataset", type=int, help="The dataset number for this simulation")
parser.add_argument("-f", "--filenr", type=int, help="The file number (within the run)")
parser.add_argument("-n", "--numevents", type=int, help="The number of events per run")
parser.add_argument("--min-energy", type=float, help="Low energy bound for MuonGun [GeV]")
parser.add_argument("--max-energy", type=float, help="Upper energy bound for MuonGun [GeV]")
parser.add_argument("--power-law-index", type=float, default=-3.0, help="Power law index (negative)")
parser.add_argument("--power-law-offset", type=float, default=150., help="Power law offset [GeV]")
parser.add_argument("--inner-cylinder", action="store_true", help="Use an inner target cylinder")
parser.add_argument("--inner-cylinder-x", type=float, default=46.29, help="Inner target cylinder origin x [m]")
parser.add_argument("--inner-cylinder-y", type=float, default=-34.88, help="Inner target cylinder origin y [m]")
parser.add_argument("--inner-cylinder-z", type=float, default=-300, help="Inner target cylinder origin z [m]")
parser.add_argument("--inner-cylinder-radius", type=float, default=150., help="Inner target cylinder radius [m]")
parser.add_argument("--inner-cylinder-length", type=float, default=500., help="Inner target cylinder length [m]")
parser.add_argument("--outer-cylinder-x", type=float, default=0., help="Outer target cylinder origin x [m]")
parser.add_argument("--outer-cylinder-y", type=float, default=0., help="Outer target cylinder origin y [m]")
parser.add_argument("--outer-cylinder-z", type=float, default=-700., help="Outer target cylinder origin z [m]")
parser.add_argument("--outer-cylinder-radius", type=float, default=1600., help="Outer target cylinder radius [m]")
parser.add_argument("--outer-cylinder-length", type=float, default=300., help="Outer target cylinder length [m]")
parser.add_argument('--kde', action='store_true', help="Use a KDE pre-scale to weight event generation towards events likely to survive to higher processing levels")
args = parser.parse_args()


print "------ settings received by a job ------ "
for k,v in sorted(vars(args).items()): 
    print( "    {0}: {1}".format(k,v))
print "------ end ------ "
#########################
# I3 
#########################

tray = I3Tray()


#########################
# FLUX 
#########################

# Flux model
model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL') #TODO Do we need a new fit to the latest CORSIKA (including SYBYLL 2.3c)???
model.flux.min_multiplicity = 1
model.flux.max_multiplicity = 1

# Check the spectral index is negative
# This is for consistency with CORSIKA scripts
# Will convert to a positve number (as MuonGun expects it) later
assert args.power_law_index <= 0., "Spectral index must be negative"
power_law_index = args.power_law_index
positive_power_law_index = -1. * power_law_index

# Muon spectrum (spectral index, break in power law, and energy range)
power_law_offset = args.power_law_offset * icecube.icetray.I3Units.GeV
min_energy = args.min_energy * icecube.icetray.I3Units.GeV
max_energy = args.max_energy * icecube.icetray.I3Units.GeV
spectrum = MuonGun.OffsetPowerLaw(
    gamma=positive_power_law_index,
    offset=power_law_offset,
    min=min_energy,
    max=max_energy,
)

#TODO Are these stored in the I frame?


#########################
# DETECTOR 
#########################

# Create an outer surface
# This will be the generation volume
outer_surface_center = icecube.dataclasses.I3Position(
    args.outer_cylinder_x,
    args.outer_cylinder_y,
    args.outer_cylinder_z,
) * icecube.icetray.I3Units.m 

outer_surface = icecube.MuonGun.Cylinder(
    length=args.outer_cylinder_length * icecube.icetray.I3Units.m,
    radius=args.outer_cylinder_radius * icecube.icetray.I3Units.m,
    center=outer_surface_center,
)

print("Outer cylinder :")
print("  Origin : %s" % outer_surface.center)
print("  Length : %s" % outer_surface.length)
print("  Radius : %s" % outer_surface.radius)


# Create the generator
# Can either use the whole detector, or set up an energy-dependent inner target surface
if args.inner_cylinder:

    inner_surface_center = icecube.dataclasses.I3Position(
        args.inner_cylinder_x,
        args.inner_cylinder_y,
        args.inner_cylinder_z,
    ) * icecube.icetray.I3Units.m

    inner_surface = icecube.MuonGun.Cylinder(
        length=args.inner_cylinder_length * icecube.icetray.I3Units.m,
        radius=args.inner_cylinder_radius * icecube.icetray.I3Units.m,
        center=inner_surface_center,
    )

    print("Inner cylinder :")
    print("  Origin : %s" % inner_surface.center)
    print("  Length : %s" % inner_surface.length)
    print("  Radius : %s" % inner_surface.radius)

    scaling = icecube.MuonGun.ConstantSurfaceScalingFunction(inner_surface)

    generator = icecube.MuonGun.EnergyDependentSurfaceInjector(
        surface=outer_surface, 
        flux=model.flux, 
        energy=spectrum,
        radius=model.radius, 
        scaling=scaling,
    )

else:

    generator = icecube.MuonGun.StaticSurfaceInjector(
        outer_surface, model.flux, spectrum, model.radius)


#TODO Are these stored in the I frame


#########################
# MUONGUN 
#########################

# Use dataset number for seeding
seed = 139005

# Set up random number generator
from globals import max_num_files_per_dataset
randomService = phys_services.I3SPRNGRandomService(
        seed = seed*2, #TODO Store to I frame
        nstreams = max_num_files_per_dataset,
        streamnum = args.filenr)
tray.context['I3RandomService'] = randomService

# Generate bundles of muons
tray.AddSegment(GenerateBundles, 'BundleGen', Generator=generator, NEvents=args.numevents, GCDFile=expandvars(args.gcdfile))

#TODO simprod-scripts/python/segments/GenerateCosmicRayMuons.py add generator module in a differnt way (see line below). Are these equivalent?
# tray.AddModule('I3MuonGun::GeneratorModule',name,Generator=num_events*generator)

# Propogate the muons in the detector
tray.AddModule('I3PropagatorModule', PropagatorServices=get_propagators(), RandomService=randomService)

# Calculate event weights
raw_weight_name = 'MuonWeight'
tray.AddModule('I3MuonGun::WeightCalculatorModule', raw_weight_name, Model=model, Generator=generator)



#########################
# APPLY KDE
#########################

weight_dict_name = "I3MCWeightDict"

if args.kde:

    import sys
    sys.path.append( expandvars('$I3_SRC/kde_filter') ) #TODO Rearrange directory to add a proper python directory structure
    from kde_filter import initialize_kde, evaluate_kde

    list_of_vars = ['en','cz','z','weights']
    kernels, kde_max_factor, kde_num_events = initialize_kde(args.numevents,list_of_vars,seed)
    tray.AddModule(evaluate_kde, 'kde', 
        Streams = [icetray.I3Frame.DAQ], 
        kernels = kernels, 
        kde_max_factor = kde_max_factor, 
        num_events = kde_num_events,
        raw_weight_name = raw_weight_name,
        weight_dict_name = weight_dict_name,
        )


#########################
# Weight dict 
#########################

# Add everything to weight dict
def write_weight_dict(frame,weight_dict_name) :

    # Create a weight dict if doesn't already exist
    # Or grab an existing one (removing it from frame ready to write the new one)
    if weight_dict_name in frame :
        weight_dict = frame[weight_dict_name]
        frame.Delete(weight_dict_name)
    else :
        from icecube.dataclasses import I3MapStringDouble
        weight_dict = I3MapStringDouble()

    # Write the event generation properties
    weight_dict["power_law_index"] = power_law_index
    weight_dict["power_law_offset"] = power_law_offset
    weight_dict["min_energy"] = min_energy
    weight_dict["max_energy"] = max_energy
    weight_dict["min_multiplicity"] = model.flux.min_multiplicity
    weight_dict["max_multiplicity"] = model.flux.max_multiplicity

    weight_dict["outer_surface_center_x"] = outer_surface.center.x
    weight_dict["outer_surface_center_y"] = outer_surface.center.y
    weight_dict["outer_surface_center_z"] = outer_surface.center.z
    weight_dict["outer_surface_length"] = outer_surface.length
    weight_dict["outer_surface_radius"] = outer_surface.radius

    weight_dict["use_inner_surface"] = int(args.inner_cylinder)
    if args.inner_cylinder :
        weight_dict["inner_surface_center_x"] = inner_surface.center.x
        weight_dict["inner_surface_center_y"] = inner_surface.center.y
        weight_dict["inner_surface_center_z"] = inner_surface.center.z
        weight_dict["inner_surface_length"] = inner_surface.length
        weight_dict["inner_surface_radius"] = inner_surface.radius

    weight_dict["use_kde"] = int(args.kde)

    # Add weight information (if KDE not run, which handles this)
    #TODO Need to better integrate this with the KDE code
    if not args.kde :
        assert frame.Has(raw_weight_name), "Muon weight '%s' not found in frame" % raw_weight_name
        raw_weight = frame[raw_weight_name].value # Weight before normalisation or KDE prescale
        num_events = float(args.numevents)
        weight = raw_weight / num_events # scale by probability of number of events in file #TODO Isn;t the generator (which is passed to the weighter) already aware of the num events?
        weight_dict["raw_weight"] = raw_weight
        weight_dict["num_events"] = num_events
        weight_dict["weight"] = weight

    # Write the weight dict
    frame.Put(weight_dict_name,weight_dict)


tray.AddModule(write_weight_dict, 'write_weight_dict', 
    Streams = [icetray.I3Frame.DAQ], 
    weight_dict_name = weight_dict_name,
)


#########################
# OUTPUT AND EXIT 
#########################

# Add an event header
tray.AddModule(
    "I3MCEventHeaderGenerator",
    "gen_header",
    Year=2009,
    DAQTime=158100000000000000,
    RunNumber=args.filenr, # Use file number so that get unique [run,event] numbers for each file in a given dataset
    EventID=1,
    IncrementEventID=True,
)

tray.AddModule('I3Writer', 'writer',
        Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation],
        filename=args.outfile)

tray.AddModule('TrashCan', 'YesWeCan')
tray.Execute()