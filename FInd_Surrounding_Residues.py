#!/sw/mcm/app/anaconda/2.3.0/envs/mdaenv/bin/python

import os
import sys
import subprocess

#subprocess.call("/hits/fast/mcm/nandekpl/software/load_anaconda.sh", shell=True)

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis
import matplotlib.cm as cmx
import matplotlib.colors as colors
#%matplotlib inline
from collections import Counter
import json


import argparse

#os.chdir("/hits/fast/mcm/nandekpl/mustafa_shared/")


parser = argparse.ArgumentParser(description='Counts the % Occurance of protein/lipid residues within a specified cutoff distance from a ligand/residue. "Command:" "python FInd_Surrounding_Residues.py -i traj_at.gro -t traj_at.xtc -fa "resname POPC" -fin "protein" -d 5 -step 10 -o occupancy.dat" ')


parser.add_argument('-i', '--input1', dest = "grofile", help='Input File name of gro file (Example: docked_ami_100ns_nowater.gro).')
parser.add_argument('-t', '--input2', dest = "trajfile", help='Input File name of gromacs trajectory file (Example: docked_ami_100ns_nowater.xtc).')
parser.add_argument('-fa', '--input3', dest = "find_around", help='Mention the selection of atoms around which you want to find residues. Specfy the selections in double quoates (" ") (Example: "resname UNK" or "protein" or "resid 458").')
parser.add_argument('-fin', '--input4', dest = "find_in", help='Mention the selection of atoms in which you want to find residues. Specfy the selections in double quoates (" ") (Example: "resname POPC" or "protein" or "resid 1-458").')
parser.add_argument('-d', '--input5', dest = "distcutoff", help='Input cutoff distance from the ligand that is taken into account in angstroms (Example: 5).')
parser.add_argument('-step', '--input6', dest = "step", help='Input iteger to specify the Interval between reading frames from MD trajectory (Example: 10). The default value is 10.', default=10)
parser.add_argument('-o', '--output1', dest = "fout", help='Name of OUTPUT DATA file', default='occupancy.dat')

#parser.add_argument('-step', action='store_true'_StoreTrueAction(option_strings=['-step'], dest='step', nargs=0, const=True, default=False, type=None, choices=None, help=None, metavar=None)

args = parser.parse_args()

PATH = os.getcwd()
print 'Current working directory is', PATH

def occurance_count(grofile, xtcfile, find_around, find_in, distcutoff, step=10):
    """Counts the occurance of protein residues within a specified cutoff from a ligand
        :Arguments:
            *grofile*
                string of the name of the GROMACS gro file
            *xtcfile*
                trajectory file
                
            *find around*
                string of the selection atoms around which you want to find residues. Specfy the selections in double quoates (" ")

            *Selection to find in*
                string of the selelction name

	     *cutoff*
                cutoff distance from the ligand that is taken into account in angstroms
             *interval*
                Input iteger to specify the Interval between reading frames from MD trajectory

        :Returns:
            *Residues and % Occupancy*
                a dictionary of residues and how many times during the simulation they are within
                a cutoff value from the ligand"""


    md_sim = MDAnalysis.Universe(grofile, xtcfile)

    frame_dict = {}
    fnum = md_sim.trajectory.n_frames
    print "Total number of frames are ", fnum
    
    skip = int(args.step)
    print "The number of frames skiped = ", skip
    firstframe_ps=None
    i = 0 
    for frame in md_sim.trajectory[::int(skip)]:
        i = i+1
        selection = md_sim.select_atoms(find_in+' and around '+distcutoff+' '+find_around)
        residue_list = [atom.residue for atom in selection]
        frame_dict[frame.time]=set(residue_list)
        if firstframe_ps == None:
            firstframe_ps = frame.time
            
    print "Number of frame read = ", i

    lastframe = frame.time
    residue_counts = Counter([item for sublist in frame_dict.values() for item in sublist])    
    #print residue_counts
    ftemp = open("temp.dat", 'w')
    #print "Residue Number Times Occupancy"
    
    for res in residue_counts:
	
        occupancy = (residue_counts[res]*100)/i
        #print res, residue_counts[res], occupancy
        
        print >> ftemp, res, occupancy
    
    ftemp.close()
    
    """
    print args.grofile
    print args.trajfile
    print args.find_around
    print args.find_in
    print args.distcutoff """
    
occurance_count (args.grofile, args.trajfile, args.find_around, args.find_in, args.distcutoff)
#occurance_count ('traj_at.gro', 'traj_at.xtc', 'protein', 'POPC', '5')

occupancy = open("occupancy.dat", 'w')

fout = open(args.fout, 'w')
print >> fout, "Residue Number %Occupancy"
os.system("sed 's/<Residue//g;s/,//g;s/>//g' temp.dat | sort -r -k3 -n > occupancy.dat")
occupancy = open("occupancy.dat", 'r')
os.system("rm temp.dat")

for line in occupancy.readlines():
    fout.write(line)

if args.fout:
    print "THE RESULTS ARE IN OUTPUT FILE:", args.fout
    #os.system("rm occupancy.dat")

 
fout.close()
occupancy.close()

