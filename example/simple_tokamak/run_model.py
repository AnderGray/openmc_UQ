#!/usr/bin/env python
import argparse
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

# XS library.
XS_LIB = inputs_dict["cross_section_path"]

# Scores to extract
scores = inputs_dict["scores"]

# File to save to
output = inputs_dict["output"]

# Generate random HDF5 file of "nuclides".
# Parallel using simple multiprocessing
N_workers = int(os.getenv('SLURM_NTASKS'))
print(f"Sampling sandy with N_workers={N_workers}")

print("Beginning NJOY processing")
with Pool(N_workers) as pool:
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
openmc_uq.run_openmc(openmc_xml_dir, random_nuc,  cross_sections_xml=XS_LIB, threads = N_workers)

# Gather results
sp = openmc.StatePoint("openmc_sim/statepoint.10.h5")

result_str=""
for score in scores:
    tally = sp.get_tally(name=score)
    df = tally.get_pandas_dataframe()
    mean = df['mean'].sum()
    std = df['std. dev.'].sum()
    tally_str = "{} {}".format(mean,std)
    result_str=results_str+tally_str

#   Store results of interest to openmc.out
with open(output, 'w') as f:
    f.write(results_str)
