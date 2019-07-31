#!/usr/bin/env python

from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units, I3Frame, OMKey
import numpy as np
import argparse, SimAnalysis

parser = argparse.ArgumentParser(description = "Find neutrino distribution from NuGen simulation")
parser.add_argument( '-n', '--minFileNum', help = "smallest file number used" )
parser.add_argument( '-N', '--maxFileNum', help = "largest file number used")
parser.add_argument( '-r', '--runType', help = "The run type (must align with GCD type)")
parser.add_argument( '-H', '--hitThresh', help = "number of hits to retain a DOM")
parser.add_argument( '-d', '--domThresh', help = "number of retained DOMs to retain a frame")
args = parser.parse_args()

# open file
infileList = []
for i in range(int(args.minFileNum), int(args.maxFileNum) + 1):
    infile = dataio.I3File('/home/dvir/workFolder/I3Files/nugen/nugenStep3/' + str(args.runType) + '/NuGen_step3_' + str(args.runType) + '_' + str(i) + '.i3.gz')
    infileList.append(infile)

successfulCombinationsList = []
 
def generateCombinations(inputNums, combination, loopStart, index, combinationFunction): 
                          
    # finished constructing combination
    if (index == len(combination)):  
        combinationFunction(combination)
        return

    i = loopStart
    while(i <= len(inputNums)-1 and len(inputNums) - i >= len(combination) - index): 
        combination[index] = inputNums[i]
        generateCombinations(inputNums, combination, i + 1, index + 1, combinationFunction) 
        i += 1

def evaluateGeometry(domHashes):
    omkeys = []
    for domHash in domHashes:
        for i in range(5):
             string = 1 + int(domHash/15) + 15*i
             domNum = int(domHash % 15)
             omkey = OMKey(string, domNum, 0)
             omkeys.append(omkey)
    
    retainedFrames, totalFrames = SimAnalysis.calculateRetainedFrames(infileList, omkeys, int(args.hitThresh), int(args.domThresh))
    
    print retainedFrames, totalFrames


evaluateGeometry(np.linspace(1,44,45))

