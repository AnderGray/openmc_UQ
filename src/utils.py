import os
from sandy import Endf6
import re
import numpy as np

from pathlib import Path

import openmc

###
#   Pattern matches nuclides with ENDF files, and returns their path
###
def get_nuclide_paths(endf_path, nuclides):
    print(f"######################################")
    print(f"##### SEARCHING FOR  ENDF FILES  #####")
    print("\n")

    file_list = np.array(os.listdir(endf_path))

    if len(file_list) == 0:
        raise Exception(f"{endf_path} is empty")

    nuclide_paths = []

    for (i, nuc) in enumerate(nuclides):
        atomic_sym = " ".join(re.findall("[a-zA-Z]+", nuc))
        mass_num = " ".join(re.findall(r'\d+', nuc))

        atomic_in = np.array([atomic_sym in files for files in file_list])
        mass_in = np.array([mass_num in files for files in file_list])

        this_file = file_list[atomic_in * mass_in]

        if len(this_file) == 0:
            nuclide_paths.append("")
            print(f"No files matched for {nuc}")
            print("\n")

        elif len(this_file) == 1:
            nuc_path = endf_path + '/' + this_file[0]
            nuclide_paths.append(nuc_path)
            print(f"Using {this_file[0]} for {nuc}")
            print("\n")

        elif len(this_file) > 1:
            shortest = min(this_file, key=len)
            nuclide_paths.append(endf_path + '/' + shortest)
            print(f"Multiple files matched for {nuc}: {this_file}")
            print(f"Using shortest: {shortest}:")
            print("\n")
            
    return nuclide_paths


###
#   Creates error and pendf files from endf file
###
def process_with_njoy(Nuclide, endf_path):

    FILE = get_nuclide_paths(endf_path, [Nuclide])[0]

    endf6 = Endf6.from_file(FILE)

    if len(endf6.mat) == 1:
        MAT = endf6.mat[0]
    else: 
        raise ValueError("ENDF file conatains multiple nuclides (MAT) numbers. Function only allows for files containing a single nuclide")

    error = endf6.get_errorr(
        err=0.001,
        irespr=0
    )['errorr33']

    pendf = endf6.get_pendf(err=0.001, verbose=True)
    
    return error, pendf


##
#   Replaces the nuclides in mat.xml with the random nuclides of NuclideStream
#   Returns openmc materials object with changed nuclides
##
def replace_nuclide_material(nuc, newNuc):
    fin = open("materials.xml", "rt")
    fout = open("materials1.xml", "wt")

    for line in fin:
	#read replace the string and write to output file
	    fout.write(line.replace(nuc, newNuc))
    
    fin.close()
    fout.close()
    os.system("mv materials1.xml materials.xml")


def replace_nuclide_tally(nuc, newNuc):
    fin = open("tallies.xml", "rt")
    fout = open("tallies1.xml", "wt")

    for line in fin:
	#read replace the string and write to output file
	    fout.write(line.replace(nuc, newNuc))
    
    fin.close()
    fout.close()
    os.system("mv tallies1.xml tallies.xml")