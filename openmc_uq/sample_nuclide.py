import os
import numpy as np
import pandas as pd
import h5py
import sandy
from sandy import Endf6
from pathlib import Path

import openmc.data

from .utils import get_nuclide_paths, get_cov


class RandomData:
    """Container for a perturbed nuclide's name and output path."""

    def __init__(self, nuclide: str, perturbed_name: str, path: Path):
        self.nuclide = nuclide
        self.perturbed_name = perturbed_name
        self.path = path


def sample_nuclide_sandy(nuclide: str, endf_path: str, seed33: int) -> RandomData:
    """
    Generate a single random nuclear data sample using SANDY's built-in sampler.
    Outputs an ACE file and converts it to HDF5 for use with OpenMC.
    """
    file = get_nuclide_paths(endf_path, [nuclide])[0]

    out_dir = Path("ND_sample").resolve()
    out_dir.mkdir(exist_ok=True)

    acetape = out_dir / f"ace_{nuclide}_rand"
    sandy_command = (
        f"python3 -m sandy.sampling {file} "
        f"--samples 1 --outname {acetape} --seed33 {seed33} "
        f"--mf 33 --acer --temperature 293.6"
    )
    os.system(sandy_command)

    data = openmc.data.IncidentNeutron.from_ace(str(acetape) + ".02c")
    data.name = f"{nuclide}_rand"

    file_out = out_dir / f"{nuclide}_rand.h5"
    data.export_to_hdf5(file_out, "w")

    return RandomData(nuclide, data.name, file_out)


def _load_svd(svd_file: Path, n_trunc: int) -> tuple[np.ndarray, np.ndarray]:
    """Load and truncate U and s from a saved SVD HDF5 file."""
    with h5py.File(svd_file, "r") as hf:
        u = np.array(hf["u"])[:, :n_trunc]
        s = np.array(hf["s"])[:n_trunc]
    return u, s


def _get_valid_mts(errorr, pendf) -> np.ndarray:
    """Return MT numbers present in both errorr and pendf, excluding MT 451."""
    mts = np.array([mt for mt in errorr.mt if mt in pendf.mt])
    return np.setdiff1d(mts, 451)  # MT 451 contains general info, not a reaction


def _apply_perturbation(xs, errorr, pendf, mat: int, mts: np.ndarray, sample: np.ndarray) -> sandy.Xs:
    """
    Build a perturbation vector from the KL sample and apply it to each MT cross section.
    Clamps resulting perturbation factors to [0, 2].
    """
    energy_grid = errorr.get_energy_grid()

    cov = get_cov(errorr)
    perturbation_factors = sample + 1.0
    perturbation_factors = np.clip(perturbation_factors, 0.0, 2.0)

    samples_df = pd.DataFrame(perturbation_factors, index=cov.data.index)
    perts = sandy.Samples(samples_df).data.T

    for mt in mts:
        p = sandy.Pert(perts[(mat, mt)].values.flatten(), index=energy_grid[1:])
        xs = xs.custom_perturbation(mat, mt, p)

    return xs.reconstruct_sums()


def sample_nuclide_KL(
    nuclide: str,
    endf_path: str,
    pre_processed_path: str,
    n_trunc: int,
    sample: np.ndarray,
) -> RandomData:
    """
    Generate a perturbed nuclear data sample using a KL (SVD-based) expansion.
    Loads pre-computed errorr/pendf/SVD files, applies the perturbation, and
    exports the result as an ACE + HDF5 file for use with OpenMC.
    """
    file = get_nuclide_paths(endf_path, [nuclide])[0]
    endf6 = Endf6.from_file(file)

    if len(endf6.mat) != 1:
        raise ValueError(
            "ENDF file contains multiple nuclides (MAT numbers). "
            "This function only supports files with a single nuclide."
        )
    mat = endf6.mat[0]

    pre_processed_path = Path(pre_processed_path).resolve()
    errorr = sandy.errorr.Errorr.from_file(pre_processed_path / f"{nuclide}.error")
    pendf  = sandy.endf6.Endf6.from_file(pre_processed_path  / f"{nuclide}.pendf")

    mts = _get_valid_mts(errorr, pendf)
    reactions = [openmc.data.reaction.REACTION_NAME[mt] for mt in mts]
    print(f"Sampling MTs: {mts}")
    print(f"Reactions:    {reactions}")

    u, s = _load_svd(pre_processed_path / f"{nuclide}_SVD.h5", n_trunc)
    kl_sample = (u @ np.sqrt(np.diag(s)) @ sample).T

    xs_perturbed = _apply_perturbation(
        sandy.Xs.from_endf6(pendf), errorr, pendf, mat, mts, kl_sample
    )

    # Export perturbed cross sections to ACE and HDF5
    out_dir = Path("sample_rand").resolve()
    out_dir.mkdir(exist_ok=True)

    pendf_tape = out_dir / f"pendf_{nuclide}_rand"
    xs_perturbed.to_endf6(pendf).to_file(pendf_tape)

    outputs = sandy.njoy.process_neutron(
        file,
        pendftape=str(pendf_tape),
        wdir=out_dir,
        keep_pendf=False,
        temperatures=[293.6],
        verbose=True,
    )

    ace_tape = out_dir / f"ace_{nuclide}_rand"
    ace_tape.write_text(outputs["ace"])
    print(f"ACE file written to '{ace_tape}'")

    data = openmc.data.IncidentNeutron.from_ace(str(ace_tape))
    data.name = f"{nuclide}-rand"

    file_out = out_dir / f"{nuclide}-rand.h5"
    data.export_to_hdf5(file_out, "w")

    return RandomData(nuclide, data.name, file_out)
