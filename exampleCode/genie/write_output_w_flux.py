from icecube import icetray, millipede, dataclasses
from I3Tray import I3Tray
from icecube.tableio import I3TableWriter, I3CSVTableService
from icecube.hdfwriter import I3HDFTableService
from math import cos

import numpy as np

import argparse

parser = argparse.ArgumentParser(description='Get Files')
parser.add_argument("-i", "--infile",default=[], type=str, nargs='+', dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_argument("-o", "--outfile",default=[], type=str, dest="OUTFILE", help="Write output to FILE")
args = parser.parse_args()

from icecube import NuFlux
flux_service  = NuFlux.makeFlux("IPhonda2014_spl_solmin")

from icecube.dataclasses import get_most_energetic_neutrino


def GetFlux(frame):

    true_neutrino = get_most_energetic_neutrino(frame["I3MCTree"])
    true_nu_energy  = true_neutrino.energy
    true_nu_coszen  = cos(true_neutrino.dir.zenith)
    norm = (frame["I3MCWeightDict"]['OneWeight'] / frame["I3MCWeightDict"]['NEvents']) * 2.

    
    if (true_neutrino.type > 0):
        nue_flux_vector  = flux_service.getFlux(dataclasses.I3Particle.NuE , true_nu_energy, true_nu_coszen) * norm*0.5/0.7
        numu_flux_vector = flux_service.getFlux(dataclasses.I3Particle.NuMu, true_nu_energy, true_nu_coszen) * norm*0.5/0.7
    else:
        nue_flux_vector  = flux_service.getFlux(dataclasses.I3Particle.NuEBar , true_nu_energy, true_nu_coszen) * norm*0.5/0.3
        numu_flux_vector = flux_service.getFlux(dataclasses.I3Particle.NuMuBar, true_nu_energy, true_nu_coszen) * norm*0.5/0.3

#    print true_nu_energy, true_nu_coszen, norm, numu_flux_vector, nue_flux_vector
    frame["I3MCWeightDict"]["no_flux"] = norm
    frame["I3MCWeightDict"]["numu_flux"] = numu_flux_vector
    frame["I3MCWeightDict"]["nue_flux"] = nue_flux_vector

    
tray = I3Tray()
tray.AddModule('I3Reader', 'reader', FilenameList = args.INFILE)

hdf = I3HDFTableService(args.OUTFILE)

tray.AddModule(GetFlux, "GetFlux", Streams = [icetray.I3Frame.DAQ])

tray.AddModule(I3TableWriter,
               tableservice = [hdf],
#               BookEverything = True,
               keys         = ['I3EventHeader','I3MCTree', 'I3MCWeightDict', 'FilterMask'], 
               SubEventStreams = ["InIceSplit"]
              )

tray.Execute()

