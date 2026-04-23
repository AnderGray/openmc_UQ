#!/usr/bin/env python
import argparse
import json
import numpy as np

from openmc_uq.utils import (
    check_covariance,
    get_nuclide_paths,
    resolve_nuclides,
)


parser = argparse.ArgumentParser()
parser.add_argument("input_file")
args = parser.parse_args()

print("Configuring from: ", args.input_file)
inputs_dict = {}
with open(args.input_file) as handle:
    inputs_dict = json.load(handle)

solver_inputs = inputs_dict["solver_inputs"]
nuclides = resolve_nuclides(
    solver_inputs.get("nuclides"),
    solver_inputs["openmc_xml_dir"],
)
nuclides = np.array(nuclides)
n_requested = len(nuclides)

endf_folder = solver_inputs["endf_dir"]
paths = np.array(get_nuclide_paths(endf_folder, nuclides))

is_not_empty = np.array([p != "" for p in paths], dtype=bool)
nuclides = nuclides[is_not_empty]
paths = paths[is_not_empty]

have_covs = np.array([check_covariance(path) for path in paths], dtype=bool)
nuclides = nuclides[have_covs]

if len(nuclides) == 0:
    raise RuntimeError(
        "No nuclides remain after filtering for available ENDF files and covariance data."
    )

if len(nuclides) < n_requested:
    print(
        f"Filtered out {n_requested - len(nuclides)} nuclides due to missing ENDF files/covariance data."
    )

inputs_dict["solver_inputs"]["nuclides"] = nuclides.tolist()

print(
    f"Writing updated solver_inputs with {len(nuclides)} filtered nuclides "
    f"to '{args.input_file}'."
)

with open(args.input_file, "w") as handle:
    json.dump(inputs_dict, handle, indent=4)
