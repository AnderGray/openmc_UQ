import os
import sandy
from sandy import Endf6
import re
import numpy as np
from scipy import linalg
from scipy.linalg import svd

from pathlib import Path
import h5py

import pandas as pd

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


###
#   Checks if the Covariance file for a cross section exists in an  endf file
###
def check_covariance(endf_file):

    isExist = os.path.exists(endf_file)

    if isExist is False:
        print(f"{endf_file} does not exist")
        return False

    tape = sandy.Endf6.from_file(endf_file)

    if 33 not in tape.mf:
        print(f"No cross section covariance info exists for {endf_file}")
        return False
    
    return True


def get_cov(error):
    
    try: 
        cov = error.get_cov()

    except:
        print("Error occured in when trying to return covariance matrix. Potentially due to negative variances. Making variances postive")

        eg = error.get_energy_grid()
        eg = pd.IntervalIndex.from_breaks(eg)  # multigroup

        # initialize global cov matrix with all MAT, MT
        ix = pd.DataFrame(error.filter_by(listmf=[31, 33]).data.keys(),
                            columns=["MAT", "MF", "MT"])[["MAT", "MT"]]
        ix["IMIN"] = ix.index * eg.size
        ix["IMAX"] = (ix.index + 1) * eg.size
        nd = ix.shape[0]
        nsize = nd * eg.size
        c = np.zeros((nsize, nsize))

        # Fill matrix
        for mat, mf, mt in error.filter_by(listmf=[31, 33]).data:
            mf33 = sandy.errorr.read_mf33(error, mat, mt)

            for mt1, cov in mf33["COVS"].items():
                ivals = ix.query("MAT==@mat & MT==@mt").squeeze()
                imin, imax = ivals.IMIN, ivals.IMAX
                jvals = ix.query("MAT==@mat & MT==@mt1").squeeze()
                jmin, jmax = jvals.IMIN, jvals.IMAX
                c[imin: imax, jmin: jmax] = cov
                if mt != mt1:
                    c[jmin: jmax, imin: imax] = cov.T

        variances = np.diag(c)
        neg_indecies = variances<0
        neg_variances = variances[variances<0]

        c[neg_indecies, neg_indecies] = np.abs(c[neg_indecies, neg_indecies])

        print("################")
        print(f"##  found {len(neg_variances)} negative variances: {neg_variances} ##")

        # Add index and columns and convert to CategoryCov
        idx = pd.MultiIndex.from_tuples(
            [(mat, mt, e) for i, (mat, mt) in ix[["MAT", "MT"]].iterrows() for e in eg],
            names=["MAT", "MT", "E"],
        )

        cov = sandy.CategoryCov(c, index=idx, columns=idx)

    return cov


def svd_expand(cov, truncation_order):
    
    u, s, vh = svd(cov, check_finite=False)

    s_ = s[:truncation_order]
    u_ = u[:, :truncation_order]

    captured_var = np.sum(s_)
    total_var = np.sum(s)

    C_test = u[:, :s_.size] @ np.diag(s_) @ u[:, :s_.size].T

    print("*******************************************************************")
    print(f" Captured {captured_var/ total_var * 100}% of the variance ")
    print("*******************************************************************")
    print('Max difference C and C_test = %.8f' % np.max(np.abs(C_test - cov)))
    print("*******************************************************************")
    print()

    return u_, s_, vh

####
#   Checks number of dimensions required to represent a nuclear data file
####
def get_number_of_dimensions(Nuclide, endf_path, total_var = 99.9):
     
    FILE = get_nuclide_paths(endf_path, [Nuclide])[0]

    endf6 = Endf6.from_file(FILE)

    if len(endf6.mat) == 1:
        MAT = endf6.mat[0]
    else: 
        raise ValueError("ENDF file conatains multiple nuclides (MAT) numbers. Function only allows for files containing a single nuclide")

    cov = get_cov(endf6.get_errorr(
        err=0.001,
        irespr=0
    )['errorr33']).data.values

    u, s, vh = svd(cov, check_finite=False)

    s_ = s[np.cumsum(s / np.sum(s)) * 100 <= total_var]
    n_t = s_.size

    print(f"Dimension reduced from {s.size} to {n_t} for {Nuclide}")
    
    return n_t, s.size

####
#   Checks number of dimensions required to represent a nuclear data file,
#   and performs initial processing (saving error, pendf, and SVD matricies)
#
####
def svd_and_save_error(Nuclide, endf_path, save_path, force_processing = False, total_var = 99.9):

    save_path = Path(save_path).resolve()
    save_path.mkdir(exist_ok=True)

    error_file = save_path/f"{Nuclide}.error"
    pendf_file = save_path/f"{Nuclide}.pendf"
    SVD_file = save_path/f"{Nuclide}_SVD.h5"

    if error_file.is_file() and pendf_file.is_file() and not force_processing:
        error = sandy.errorr.Errorr.from_file(error_file)
        pendf = sandy.endf6.Endf6.from_file(pendf_file)

    else:
        error, pendf = process_with_njoy(Nuclide, endf_path)
        error.to_file(error_file)
        pendf.to_file(pendf_file)


    if SVD_file.is_file() and not force_processing:

        hf = h5py.File(SVD_file, 'r')
        u = np.array(hf['u'])
        s = np.array(hf['s'])
        vh = np.array(hf['vh'])

        hf.close()

    else: 
        print("****** Performing SVD ********")
        cov = get_cov(error).data.values
        u, s, vh = svd(cov, check_finite=False)

        hf = h5py.File(SVD_file, 'w')

        hf.create_dataset('u', data=u)
        hf.create_dataset('s', data=s)
        hf.create_dataset('vh', data=vh)

        hf.close()

    s_ = s[np.cumsum(s / np.sum(s)) * 100 < total_var]
    n_t = s_.size

    MTs = error.mt                 # Get available MT numbers
    MTs_pendf = pendf.mt

    is_mt_in_pendf = [mt in MTs_pendf for mt in MTs]

    MTs = np.array(MTs)[is_mt_in_pendf]
    
    #MTs.remove(451)                 
    MTs = np.setdiff1d(MTs, 451)        # Recentlty added MT, contains general information

    REACTION = openmc.data.reaction.REACTION_NAME
    reactions = [REACTION[mt] for mt in MTs]

    print("###########################")
    print("#                         #")
    print(f"Sampling mts: {MTs}")
    print(f"Reactions: {reactions}")
    print("#                         #")
    print("###########################")

    print(f"Dimension reduced from {s.size} to {n_t} for {Nuclide}")
    
    return n_t, s.size