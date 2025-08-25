#!/usr/bin/env python
import argparse
from glob import glob
import json
from multiprocessing import Pool
import os
from pathlib import Path
import sys
import openmc
import openmc_uq

parser = argparse.ArgumentParser()
parser.add_argument("input_file")
args = parser.parse_args()

print("Configuring from: ", args.input_file)
inputs_dict={}
with open(args.input_file) as handle:
    inputs_dict = json.load(handle)

# Random seeds (one for each nuclide)
seeds = inputs_dict["seeds"]

# Nuclides to sample
nuclides = inputs_dict["nuclides"]

# ENDF directory. Used to generate random data.
endf_path = inputs_dict["endf_dir"]

# Directory of openmc model.
openmc_xml_dir = inputs_dict["openmc_xml_dir"]

# Local directory name to run in
run_dir = inputs_dict["run_dir"]

# XS library.
XS_LIB = inputs_dict["cross_section_path"]

# Scores to extract
scores = inputs_dict["scores"]

# File to save to
output = inputs_dict["output"]

# Number threads / workers / cores
n_cores = int(inputs_dict["n_cores"])
n_mpi = int(inputs_dict["n_mpi"])

# Generate random HDF5 file of "nuclides".
# Parallel using simple multiprocessing
print(f"Sampling sandy with n_cores={n_cores}")

print("Beginning NJOY processing")
with Pool(n_cores) as pool:
    random_nuc = []
    for (i, nuc) in enumerate(nuclides):

        func_args = (nuc, endf_path, int(seeds[i]))
        RN = pool.apply_async(openmc_uq.sample_nuclide_sandy, func_args)
        random_nuc.append(RN)

    for r in random_nuc:
        r.wait()

    for (i, r) in enumerate(random_nuc):
        random_nuc[i] = r.get()

#  Run openmc with random files
openmc_uq.run_openmc(openmc_xml_dir, random_nuc, cross_sections_xml=XS_LIB, threads=n_cores, n_mpi = n_mpi,run_dir=run_dir)

# Lift tallies from statepoint file
n_max_batches=0
stem=run_dir+'/statepoint.'
template=stem+'*.h5'
files=glob(template)
file_to_open=""
for filename in files:
    root=filename.replace('.h5','')
    n_batches=int(root.replace(stem,''))
    if n_batches > n_max_batches:
        n_max_batches=n_batches
        file_to_open=filename

print("Opening statepoint file: ",file_to_open)
sp = openmc.StatePoint(file_to_open)

result_str=""
for score in scores:
    tally = sp.get_tally(name=score)
    df = tally.get_pandas_dataframe()
    mean = df['mean'].sum()
    std = df['std. dev.'].sum()
    tally_str = "{} {} ".format(mean,std)
    result_str=result_str+tally_str

#   Store results of interest to openmc.out
with open(output, 'w') as f:
    f.write(result_str)
