#!/usr/bin/env python

import h5py
import numpy as np
import math
import argparse
import os, sys

parser = argparse.ArgumentParser(description='Get Files')
parser.add_argument("-i", "--infile",default=[], type=str, nargs='+', dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")

args = parser.parse_args()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

day_sec = 60.*60.*24.
yr_sec = day_sec*365. 

plotting_dict = { "energy":  {"log": True,
                              "bounds": [0.0,4.0],
                              "label": "Log(E) (GeV)",
                              "title": "Energy"},
                  "zenith":  {"log": False,
                              "bounds": [-1.0,1.0],
                              "label": "cos(zenith)",
                              "title": "Zenith"},
                  "azimuth": {"log": False,
                              "bounds": [0.0, math.pi],
                              "label": "azimuth",
                              "title": "azimuth"},
                  "x":       {"log": False,
                              "bounds": [-300.,300.],
                              "label": "x (m)",
                              "title": "X"},
                  "y":       {"log": False,
                              "bounds": [-300.,300],
                              "label": "y (m)",
                              "title": "Y"},
                  "z":       {"log": False,
                              "bounds": [-500, 0.],
                              "label": "z (m)",
                              "title": "Z"},

    }

for f in args.INFILE:
    flv, pdg, runid, blah = os.path.basename(f).split("_")

    mc_dict=dict()
    try: 
        mc = h5py.File(f, 'r')
    except:
        print >>sys.stderr, "File %s is corrupt!" % (f)
        continue 

    if "I3MCTree" not in mc.keys():
        print >>sys.stderr, "File %s is empty!" % (f)
        continue
            
    #Get the neutrino's only 
    primary_bool = mc["I3MCTree"]["shape"] == 10
    #We only want events that pass deepcore filter
    filter_bool=mc["FilterMask"]["DeepCoreFilter_13"][:,0] == True

    #get useful vars that passed above cuts
    mc_dict["energy"]=mc["I3MCTree"]["energy"][primary_bool][filter_bool]
    mc_dict["zenith"]=np.cos(mc["I3MCTree"]["zenith"][primary_bool][filter_bool])
    mc_dict["azimuth"]=mc["I3MCTree"]["azimuth"][primary_bool][filter_bool]
    mc_dict["x"]=mc["I3MCTree"]["x"][primary_bool][filter_bool]
    mc_dict["y"]=mc["I3MCTree"]["y"][primary_bool][filter_bool]
    mc_dict["z"]=mc["I3MCTree"]["z"][primary_bool][filter_bool]    

    weight=mc["I3MCWeightDict"]["nue_flux"][filter_bool]
    
    #what is the total livetime of the entire file?
    w=sum(weight[:])
    w2=sum(weight[:]**2)
    livetime=w/w2
    print "File: %s Weights: %f weights^2: %f livetime: %f" % (f, w, w2, livetime/day_sec)

    #onto plotting:

    # #energy/bin
    # bins = 10.**np.linspace(0.0, 4.0, 81)
    # hist_w, b = np.histogram(energy, bins=bins, weights=weight)
    # hist_w2, b = np.histogram(energy, bins=bins, weights=weight**2)

    # plt.figure()
    # plt.step(bins, np.append([0], hist_w/hist_w2), color="blue")
    # plt.xscale("log")
    # plt.xlabel("Log(E)")
    # plt.title(flv+" "+runid+" Energy Livetime/bin")
    # plt.savefig(flv+"_"+pdg+"_"+runid+"_energy_livetime.png")
    # plt.close()

    # #zenith/bin
    # bins = np.linspace(-1.0, 1.0, 81)
    # hist_w, b = np.histogram(np.cos(zenith), bins=bins, weights=weight)
    # hist_w2, b = np.histogram(np.cos(zenith), bins=bins, weights=weight**2)

    # plt.figure()
    # plt.step(bins, np.append([0], hist_w/hist_w2), color="blue")
    # plt.xscale("log")
    # plt.xlabel("cos(zenith)")
    # plt.title(flv+" "+runid+" Zenith Livetime/bin")
    # plt.savefig(flv+"_"+pdg+"_"+runid+"_zenith_livetime.png")
    # plt.close()

    # #azimuth/bin
    # bins = np.linspace(0.0, math.pi, 81)
    # hist_w, b = np.histogram(azimuth, bins=bins, weights=weight)
    # hist_w2, b = np.histogram(azimuth, bins=bins, weights=weight**2)

    # plt.figure()
    # plt.step(bins, np.append([0], hist_w/hist_w2), color="blue")
    # plt.xscale("log")
    # plt.xlabel("azimuth")
    # plt.title(flv+" "+runid+" Azimuth Livetime/bin")
    # plt.savefig(flv+"_"+pdg+"_"+runid+"_azimuth_livetime.png")
    # plt.close()

    for var in mc_dict:
        if plotting_dict[var]["log"]:
            bins = 10.**np.linspace(plotting_dict[var]["bounds"][0], plotting_dict[var]["bounds"][1], 81)
        else:
            bins = np.linspace(plotting_dict[var]["bounds"][0], plotting_dict[var]["bounds"][1], 81)

        hist_w, b = np.histogram(mc_dict[var], bins=bins, weights=weight)
        hist_w2, b = np.histogram(mc_dict[var], bins=bins, weights=weight**2)
        
        plt.figure()
        plt.step(bins, np.append([0], hist_w/hist_w2), color="blue")
        if plotting_dict[var]["log"]:
            plt.xscale("log")     
        plt.xlabel(plotting_dict[var]["label"])
        plt.title(flv+" "+runid+" "+" "+plotting_dict[var]["title"]+" Livetime/bin")
        plt.savefig(flv+"_"+pdg+"_"+runid+"_"+plotting_dict[var]["title"]+"_livetime.png")
        plt.close()

        plt.figure()
        plt.hist(mc_dict[var], bins=bins, weights=weight, color="blue")
        if plotting_dict[var]["log"]:
            plt.xscale("log")     
        plt.xlabel(plotting_dict[var]["label"])
        plt.title(flv+" "+runid+" "+" "+plotting_dict[var]["title"])
        plt.savefig(flv+"_"+pdg+"_"+runid+"_"+plotting_dict[var]["title"]+".png")
        plt.close()
    
    mc.close()
