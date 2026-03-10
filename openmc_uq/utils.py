import os
import re
import h5py
import numpy as np
import pandas as pd
import sandy
from sandy import Endf6
import openmc.data.reaction
from pathlib import Path
from scipy.linalg import svd
import xml.etree.ElementTree as ET


def _unique_preserve_order(values: list[str]) -> list[str]:
    """Return unique values in first-seen order."""
    unique = []
    seen = set()
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        unique.append(value)
    return unique


def get_nuclides_from_materials_xml(materials_xml: str | Path) -> list[str]:
    """
    Parse `materials.xml` and return unique nuclide names in first-seen order.
    """
    materials_xml = Path(materials_xml).resolve()
    if not materials_xml.is_file():
        raise FileNotFoundError(
            f"Cannot auto-populate nuclides: '{materials_xml}' was not found."
        )

    try:
        root = ET.parse(materials_xml).getroot()
    except ET.ParseError as exc:
        raise ValueError(f"Failed to parse '{materials_xml}': {exc}") from exc

    nuclides = []
    for nuclide in root.findall(".//nuclide"):
        name = nuclide.get("name")
        if name:
            nuclides.append(name)

    nuclides = _unique_preserve_order(nuclides)
    if len(nuclides) == 0:
        raise ValueError(
            f"Cannot auto-populate nuclides: no '<nuclide ... name=\"...\"/>' entries "
            f"were found in '{materials_xml}'."
        )

    return nuclides


def resolve_nuclides(nuclides, openmc_xml_dir: str | Path) -> list[str]:
    """
    Resolve nuclides from explicit input or fallback to openmc_xml_dir/materials.xml.

    Fallback is triggered when `nuclides` is missing/null/empty list.
    """
    if nuclides is None or (isinstance(nuclides, list) and len(nuclides) == 0):
        materials_xml = Path(openmc_xml_dir).resolve() / "materials.xml"
        resolved = get_nuclides_from_materials_xml(materials_xml)
        print(
            f"No nuclides were specified. Loaded {len(resolved)} nuclides from "
            f"'{materials_xml}'."
        )
        return resolved

    if not isinstance(nuclides, list):
        raise TypeError(
            "Invalid `nuclides` input: expected a list, null, or missing key."
        )

    resolved = []
    for i, nuclide in enumerate(nuclides):
        if not isinstance(nuclide, str) or len(nuclide.strip()) == 0:
            raise ValueError(
                f"Invalid nuclide at index {i}: expected a non-empty string, got {nuclide!r}."
            )
        resolved.append(nuclide)

    return resolved

def get_nuclide_paths(endf_path: str, nuclides: list[str]) -> list[str]:
    """
    Pattern-match nuclide names against files in endf_path and return their full paths.
    Returns an empty string for any nuclide with no match.
    """
    print("Searching for ENDF files...")

    file_list = np.array(os.listdir(endf_path))
    if len(file_list) == 0:
        raise FileNotFoundError(f"{endf_path} is empty")

    nuclide_paths = []

    for nuclide in nuclides:
        atomic_sym = "".join(re.findall("[a-zA-Z]+", nuclide))
        mass_num   = "".join(re.findall(r"\d+", nuclide))

        matches = file_list[
            np.array([atomic_sym in f for f in file_list]) &
            np.array([mass_num   in f for f in file_list])
        ]

        if len(matches) == 0:
            print(f"  No files matched for {nuclide}")
            nuclide_paths.append("")

        elif len(matches) == 1:
            path = str(Path(endf_path) / matches[0])
            print(f"  Using {matches[0]} for {nuclide}")
            nuclide_paths.append(path)

        else:
            shortest = min(matches, key=len)
            print(f"  Multiple files matched for {nuclide}: {matches}. Using shortest: {shortest}")
            nuclide_paths.append(str(Path(endf_path) / shortest))

    return nuclide_paths


def process_with_njoy(nuclide: str, endf_path: str):
    """
    Process an ENDF file with NJOY to produce errorr and pendf outputs.
    The ENDF file must contain exactly one nuclide (MAT number).
    """
    file = get_nuclide_paths(endf_path, [nuclide])[0]
    endf6 = Endf6.from_file(file)

    if len(endf6.mat) != 1:
        raise ValueError(
            "ENDF file contains multiple nuclides (MAT numbers). "
            "This function only supports files with a single nuclide."
        )

    error = endf6.get_errorr(err=0.001, irespr=0)["errorr33"]
    pendf = endf6.get_pendf(err=0.001, verbose=True)

    return error, pendf


def _replace_in_xml(filepath: str, old: str, new: str) -> None:
    """Replace all occurrences of `old` with `new` in an XML file, in place."""
    path = Path(filepath)
    path.write_text(path.read_text().replace(old, new))


def replace_nuclide_material(old_nuclide: str, new_nuclide: str) -> None:
    """Replace a nuclide name in materials.xml."""
    _replace_in_xml("materials.xml", old_nuclide, new_nuclide)


def replace_nuclide_tally(old_nuclide: str, new_nuclide: str) -> None:
    """Replace a nuclide name in tallies.xml."""
    _replace_in_xml("tallies.xml", old_nuclide, new_nuclide)



def check_covariance(endf_file: str) -> bool:
    """Check if covariance data for a cross section exists in an ENDF file."""
    if not os.path.exists(endf_file):
        print(f"{endf_file} does not exist")
        return False

    tape = sandy.Endf6.from_file(endf_file)

    if 33 not in tape.mf:
        print(f"No cross section covariance info exists for {endf_file}")
        return False

    return True


def _build_cov_matrix(error) -> np.ndarray:
    """Build the full covariance matrix from an errorr object."""
    eg = pd.IntervalIndex.from_breaks(error.get_energy_grid())

    ix = pd.DataFrame(
        error.filter_by(listmf=[31, 33]).data.keys(),
        columns=["MAT", "MF", "MT"]
    )[["MAT", "MT"]]
    ix["IMIN"] = ix.index * eg.size
    ix["IMAX"] = (ix.index + 1) * eg.size

    nsize = ix.shape[0] * eg.size
    c = np.zeros((nsize, nsize))

    for mat, mf, mt in error.filter_by(listmf=[31, 33]).data:
        mf33 = sandy.errorr.read_mf33(error, mat, mt)
        for mt1, cov in mf33["COVS"].items():
            ivals = ix.query("MAT==@mat & MT==@mt").squeeze()
            jvals = ix.query("MAT==@mat & MT==@mt1").squeeze()
            c[ivals.IMIN:ivals.IMAX, jvals.IMIN:jvals.IMAX] = cov
            if mt != mt1:
                c[jvals.IMIN:jvals.IMAX, ivals.IMIN:ivals.IMAX] = cov.T

    return c, ix, eg, mat


def _fix_negative_variances(c: np.ndarray) -> np.ndarray:
    """Fix negative variances on the diagonal by taking their absolute value."""
    variances = np.diag(c)
    neg_indices = variances < 0
    neg_variances = variances[neg_indices]

    c[neg_indices, neg_indices] = np.abs(c[neg_indices, neg_indices])

    print(f"Found {len(neg_variances)} negative variances: {neg_variances}")
    return c


def get_cov(error) -> sandy.CategoryCov:
    """Return the covariance matrix, correcting negative variances if needed."""
    try:
        return error.get_cov()
    except Exception:
        print(
            "Error retrieving covariance matrix (possibly negative variances). "
            "Falling back to manual construction with absolute variances."
        )

    c, ix, eg, mat = _build_cov_matrix(error)
    c = _fix_negative_variances(c)

    idx = pd.MultiIndex.from_tuples(
        [(mat, mt, e) for _, (mat, mt) in ix[["MAT", "MT"]].iterrows() for e in eg],
        names=["MAT", "MT", "E"],
    )

    return sandy.CategoryCov(c, index=idx, columns=idx)


def svd_expand(cov: np.ndarray, truncation_order: int):
    """Perform SVD and truncate to the given order, reporting captured variance."""
    u, s, vh = svd(cov, check_finite=False)

    s_truncated = s[:truncation_order]
    u_truncated = u[:, :truncation_order]

    captured_pct = np.sum(s_truncated) / np.sum(s) * 100
    c_reconstructed = u_truncated @ np.diag(s_truncated) @ u_truncated.T
    max_diff = np.max(np.abs(c_reconstructed - cov))

    print(f"Captured {captured_pct:.4f}% of the variance")
    print(f"Max reconstruction error: {max_diff:.8f}")

    return u_truncated, s_truncated, vh


def get_number_of_dimensions(nuclide: str, endf_path: str, total_var: float = 99.9) -> tuple[int, int]:
    """Return the number of SVD dimensions needed to capture `total_var`% of variance."""
    file = get_nuclide_paths(endf_path, [nuclide])[0]
    endf6 = Endf6.from_file(file)

    if len(endf6.mat) != 1:
        raise ValueError(
            "ENDF file contains multiple nuclides (MAT numbers). "
            "This function only supports files with a single nuclide."
        )

    cov = get_cov(endf6.get_errorr(err=0.001, irespr=0)["errorr33"]).data.values
    u, s, vh = svd(cov, check_finite=False)

    n_truncated = np.sum(np.cumsum(s / np.sum(s)) * 100 <= total_var)
    print(f"Dimension reduced from {s.size} to {n_truncated} for {nuclide}")

    return n_truncated, s.size


def svd_and_save_error(
    nuclide: str,
    endf_path: str,
    save_path: str,
    force_processing: bool = False,
    total_var: float = 99.9,
) -> tuple[int, int]:
    """
    Process a nuclide's ENDF file: compute or load errorr/pendf outputs and SVD,
    then return the truncated and full SVD dimension counts.
    """
    save_path = Path(save_path).resolve()
    save_path.mkdir(exist_ok=True)

    error_file = save_path / f"{nuclide}.error"
    pendf_file = save_path / f"{nuclide}.pendf"
    svd_file   = save_path / f"{nuclide}_SVD.h5"

    # Load or compute errorr/pendf
    if error_file.is_file() and pendf_file.is_file() and not force_processing:
        error = sandy.errorr.Errorr.from_file(error_file)
        pendf = sandy.endf6.Endf6.from_file(pendf_file)
    else:
        error, pendf = process_with_njoy(nuclide, endf_path)
        error.to_file(error_file)
        pendf.to_file(pendf_file)

    # Load or compute SVD
    if svd_file.is_file() and not force_processing:
        with h5py.File(svd_file, "r") as hf:
            u  = np.array(hf["u"])
            s  = np.array(hf["s"])
            vh = np.array(hf["vh"])
    else:
        print("Performing SVD...")
        cov = get_cov(error).data.values
        u, s, vh = svd(cov, check_finite=False)

        with h5py.File(svd_file, "w") as hf:
            hf.create_dataset("u",  data=u)
            hf.create_dataset("s",  data=s)
            hf.create_dataset("vh", data=vh)

    n_truncated = np.sum(np.cumsum(s / np.sum(s)) * 100 < total_var)

    # Determine valid MTs present in both error and pendf
    mts = np.array([mt for mt in error.mt if mt in pendf.mt])
    mts = np.setdiff1d(mts, 451)  # MT 451 contains general info, not a reaction

    reaction_names = [openmc.data.reaction.REACTION_NAME[mt] for mt in mts]
    print(f"Sampling MTs: {mts}")
    print(f"Reactions:    {reaction_names}")
    print(f"Dimension reduced from {s.size} to {n_truncated} for {nuclide}")

    return n_truncated, s.size
