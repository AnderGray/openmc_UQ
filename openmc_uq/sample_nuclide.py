# %% Imports
import os
import openmc.data
from pathlib import Path
from .utils import get_nuclide_paths, get_cov
import sandy
from sandy import Endf6
import numpy as np

import pandas as pd
import h5py

class RandomData:
    def __init__(self, nuclide, perturbed_name, path):
        
        self.nuclide = nuclide
        self.perturbed_name = perturbed_name
        self.path = path


def sample_nuclide_sandy(Nuclide, endf_path, seed33):
    
    FILE = get_nuclide_paths(endf_path, [Nuclide])[0]


    out_dir = Path(f"ND_sample").resolve()
    out_dir.mkdir(exist_ok=True)

    acetape = os.path.join(out_dir, f"ace_{Nuclide}_rand")

    sandy_command = f"python3 -m sandy.sampling {FILE} --samples 1 --outname {acetape} --S33 {seed33} --mf 33 --acer --temperature 293.6"
    os.system(sandy_command)

    fileOut = out_dir / f"{Nuclide}_rand.h5"

    data = openmc.data.IncidentNeutron.from_ace(acetape + ".02c")
    data.name = f"{Nuclide}_rand"
    data.export_to_hdf5(fileOut, "w")

    random_nuc = RandomData(Nuclide, data.name, fileOut)

    return random_nuc

def sample_nuclide_KL(Nuclide, endf_path, pre_processed_path, n_trunc, sample):

    FILE = get_nuclide_paths(endf_path, [Nuclide])[0]
    endf6 = Endf6.from_file(FILE)

    if len(endf6.mat) == 1:
        MAT = endf6.mat[0]
    else: 
        raise ValueError("ENDF file conatains multiple nuclides (MAT) numbers. Function only allows for files containing a single nuclide")
    
    pre_processed_path = Path(pre_processed_path).resolve()
    error_file = pre_processed_path/f"{Nuclide}.error"
    pendf_file = pre_processed_path/f"{Nuclide}.pendf"
    SVD_file = pre_processed_path/f"{Nuclide}_SVD.h5"

    errorr = sandy.errorr.Errorr.from_file(error_file)
    pendf = sandy.endf6.Endf6.from_file(pendf_file)
    
    Es = errorr.get_energy_grid()   # Get energy grid 
    MTs = errorr.mt                 # Get available MT numbers
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

    hf = h5py.File(SVD_file, 'r')
    u = np.array(hf['u'])
    s = np.array(hf['s'])
    # vh = np.array(hf['vh'])

    hf.close()

    s_ = s[:n_trunc]
    u_ = u[:, :n_trunc]
    #sample = np.random.rand(n_trunc, 1)
    samples_SN = u_ @ np.sqrt(np.diag(s_)) @ sample      # Should be  u_ @ np.sqrt(np.diag(s_)) @ sample ??
    samples_SN = np.array(samples_SN).T

    samples_SN = samples_SN + 1.0

    # % Insert perturbation
    xs_new = sandy.Xs.from_endf6(pendf)


    # Insert perturbations
    print()
    print(" *************** Perturbing XS ***************")
    print()

    # samples_SN += 1
    lower_bound = samples_SN > 0
    upper_bound = samples_SN < 2
    samples_SN = np.where(lower_bound, samples_SN, 0)
    samples_SN = np.where(upper_bound, samples_SN, 2)

    cov = get_cov(errorr)

    index_data = cov.data.index
    # samples = pd.DataFrame(samples_SN, index=index, columns=columns)
    samples = pd.DataFrame(samples_SN, index=index_data)
    perts = sandy.Samples(samples).data.T

    for mt in MTs:

        p = sandy.Pert(perts[(MAT, mt)].values.flatten(), index = Es[1:])
        xs_new = xs_new.custom_perturbation(MAT, mt, p)
    
    xspert = xs_new.reconstruct_sums()
    
    # Export to ACE
    print(f"******* Generating ACE and HDF5 for {Nuclide} ********")
    out_dir = Path(f"sample_rand").resolve()
    out_dir.mkdir(exist_ok=True)

    pendftape = os.path.join(out_dir, f"pendf_{Nuclide}_rand")
    xspert.to_endf6(pendf).to_file(pendftape)
    outputs = sandy.njoy.process_neutron(
                FILE,
                pendftape=pendftape,
                wdir=out_dir,
                keep_pendf=False,
                temperatures=[293.6],
                verbose=True,
                )

    acetape = os.path.join(out_dir, f"ace_{Nuclide}_rand")

    with open(acetape, "w") as f:    
        print(f"writing to file '{f}'")
        f.write(outputs["ace"])

    fileOut = out_dir / f"{Nuclide}-rand.h5"

    data = openmc.data.IncidentNeutron.from_ace(acetape)
    data.name = f"{Nuclide}-rand"
    data.export_to_hdf5(fileOut, "w")

    random_nuc = RandomData(Nuclide, data.name, fileOut)
    
    return random_nuc