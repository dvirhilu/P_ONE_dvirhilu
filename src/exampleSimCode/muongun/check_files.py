#!/usr/bin/env python

# This file has very small modifications from the original example script. It adds MMC for NuMu files.

from optparse import OptionParser
from os.path import expandvars
import os, random

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-l", "--level",default=None,
                  dest="LEVEL", help="Level to check")
parser.add_option("-b" , "--base-folder", dest = "BASEFOLDER", 
                default = "/home/terliuk/projects/rpp-dgrant/terliuk/simulations/muongun/", 
                help="Base folder for simulations")
parser.add_option("-r", "--run", default=None,dest="RUN",
                help = "Run number" )
parser.add_option("-t", "--tmpdir", default="/scratch/terliuk/sim/muongun/",dest="TMPDIR",
                help = "" )
(options,args) = parser.parse_args()


print "------ settings received by a job ------ "
for k,v in sorted(vars(options).items()): 
    print( "    {0}: {1}".format(k,v))
print "------ end ------ "

check_folder = os.path.join(options.BASEFOLDER, str(options.LEVEL), str(options.RUN))
print("Folder: ", check_folder)
allfnames = os.listdir(check_folder)
allfnames.sort()
print "Total files found : ", len(allfnames)
last_index = allfnames[-1].split(".i3")[0].split(".")[-1]
print "last index : ", last_index
print "missing files " , int(last_index)+1 - len(allfnames)

file_prefix = allfnames[-1].split(".i3")[0].split(last_index)[0]
file_extension = allfnames[-1].split(".i3")[1].split(last_index)[0]
print "File format : ", file_prefix, file_extension
print "example filename : ", "%s%0.6i.i3%s"%(file_prefix, 0, file_extension)

file_params = { "fname" : [], 
                "size" : []}
for i in xrange(0, int(last_index)+1):
    fname = "%s%0.6i.i3%s"%(file_prefix, i, file_extension)
    file_params['fname'].append(fname)
    fullpath = os.path.join(check_folder, fname) 
    if os.path.isfile(fullpath):
        file_params['size'].append(os.path.getsize(fullpath) )
    else : 
        file_params['size'].append(-1)
tot_fsize = 0.0
tot_outfiles = 0
for s in file_params['size']:
    if s!=-1:
        tot_fsize+=s
        tot_outfiles+=1
avg_fsize = tot_fsize/tot_outfiles
sq_fsize_diff = 0.0
for s in file_params['size']:
    if s!=-1: 
        sq_fsize_diff+=(s - avg_fsize)**2
stddev = ( sq_fsize_diff/(tot_outfiles  -1) )**0.5
print "Average size : %0.3f += %0.3f MB (%i files)"%(1.0*avg_fsize/1024/1024, 1.*stddev/1024/1024, tot_outfiles)
small_file_ind = []
for i in xrange(0, int(last_index)+1):
    if file_params['size'][i]!=-1:
        if (avg_fsize - file_params['size'][i]) > 5.0*stddev:
            print "File: %s is too small : %0.3f MB" % (file_params['fname'][i], 
                                                     1.*file_params['size'][i]/1024/1024 ) 
            small_file_ind.append(i)
print "Total small files : ", len(small_file_ind)
missing_file_ind = []
for i in xrange(0,int(last_index)+1):
    if file_params['size'][i]==-1:
        if os.path.isfile(os.path.join(options.TMPDIR, file_params['fname'][i]) ):
            print "File found in temporary directory: ",file_params['fname'][i]
        else: 
            print "File missing completely : ", file_params['fname'][i]
            missing_file_ind.append(i)
script_name_dict = {"step1" : "step1", 
                    "step2" : "step2", 
                    "Level2" : "step3"}
if len(small_file_ind) > 0:
    array_string=",".join([str(i+1) for i in small_file_ind])    
    resubmit_string="sbatch --array=%s MG_%s_submit.sh %i"%(array_string, 
                                                            script_name_dict[options.LEVEL], 
                                                            int(options.RUN) )
    print "\t\t\tTo resubmit small files do :"
    print resubmit_string   
else:
    print "No small files found"
if len(missing_file_ind)>0:
    array_string=",".join([str(i+1) for i in missing_file_ind])
    resubmit_string="sbatch --array=%s MG_%s_submit.sh %i"%(array_string, 
                                                            script_name_dict[options.LEVEL], 
                                                            int(options.RUN) )

    print "\t\t\tTo resubmit missing files do :"
    print resubmit_string
else:
    print "All files were found in output or temporary directory"


                 



