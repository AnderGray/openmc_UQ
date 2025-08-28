#!/usr/bin/env python
import argparse
import json
import os
import numpy as np

from multiprocessing import Pool

from openmc_uq.utils import get_nuclide_paths, check_covariance, svd_and_save_error

parser = argparse.ArgumentParser()
parser.add_argument("input_file")
args = parser.parse_args()

print("Configuring from: ", args.input_file)
inputs_dict={}
with open(args.input_file) as handle:
    inputs_dict = json.load(handle)

nuclides =  inputs_dict["solver_inputs"]["nuclides"]

nuclides = np.array(nuclides)

endf_folder = inputs_dict["solver_inputs"]["endf_dir"]
save_folder = inputs_dict["solver_inputs"]["pre_processed_dir"]

paths = get_nuclide_paths(endf_folder, nuclides)
paths = np.array(paths)

is_not_empty = np.array([not p =='' for p in paths])
nuclides = nuclides[is_not_empty]
paths = paths[is_not_empty]

have_covs = [check_covariance(file) for file in paths]

nuclides = nuclides[have_covs]
paths = paths[have_covs]

N_workers = int(os.getenv('SLURM_NTASKS'))
print(f"Performing svd with N_workers={N_workers}")

print("Beginning NJOY processing")

with Pool(N_workers) as pool:
    async_results = []
    for (i, nuc) in enumerate(nuclides):

        func_args = (nuc, endf_folder, save_folder, False)
        RN = pool.apply_async(svd_and_save_error, func_args)
        async_results.append(RN)

    for r in async_results:
        r.wait()

    results = [r.get() for r in async_results]  # list of tuples

reduced_dims, full_dims = map(list, zip(*results))


print()
print()
print("---------------------------------------------------------------------")
print(f"Nuclides {nuclides} have been sucessfully processed and have the following reduced dimensions")
print(f"reduced_dims: {reduced_dims}")
print(f"full_dims: {full_dims}")
print("---------------------------------------------------------------------")
print()
print()
print("########## TODO ADD ADDITIONAL WARNING THAT SOME NUCLIDES MAY HAVE BEEN REMOVED FROM THE LIST IF NO COVARIANCE FOUND ############")
print()
print()

# # --- save results into JSON ---
# inputs_dict["nuclide_dimensions"] = {
#     "nuclides": list(nuclides),
#     "reduced_dims_nuclides": reduced_dims,
#     "full_dims_nuclides": full_dims
# }

if inputs_dict["dimension_reduction"]:
    nuclide_dim_sim = reduced_dims
else:
    nuclide_dim_sim = full_dims

inputs_dict["solver_inputs"]["dimensions"] = {
    "dims": nuclide_dim_sim,
    }

# Save back to file
with open(args.input_file, "w") as f:
    json.dump(inputs_dict, f, indent=4)  # indent=4 makes it pretty
