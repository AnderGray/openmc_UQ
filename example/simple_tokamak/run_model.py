import os
from pathlib import Path
import sys
from multiprocessing import Pool
import openmc


top_dir = "/home/ir-gray4/work/sims/openmc_UQ"        # Top directory of the python scripts
sys.path.insert(1, f'{top_dir}/src')
from sample_nuclide import sample_nuclide_sandy
from run_openmc import run_openmc


###
#   
###
nuclides =  ["Pb204", "Pb206", "Pb207", "Pb208", "Li6", "Li7", "W180", "W182", "W183", "W184"]  # Nuclide you wish to sample

endf_folder = f"/home/ir-gray4/work/NuclearData/tendl_2019_endf/neutron"        # ENDF directory. Used to generate random data. Change.
endf_path = str(Path(endf_folder).resolve())

openmc_problem_dir = f"{top_dir}/example/openmc_geometries/Simple_tokamak"      # Directory of openmc model. Change.
openmc_xml_dir = str(Path(openmc_problem_dir).resolve())

XS_LIB = "/home/ir-gray4/work/NuclearData/tendl-2019-hdf5/cross_sections.xml"   # XS library. Change. Make sure same as selected ENDF library


###
#   Julia script enters random seeds here. Needs to be same number as number of nuclides
###
seed = [{{{   :X1   }}}, {{{   :X2   }}}, {{{   :X3   }}}, {{{   :X4   }}}, {{{   :X5   }}}, {{{   :X6   }}}, {{{   :X7   }}}, {{{   :X8   }}}, {{{   :X9   }}}, {{{   :X10   }}}, ]


###
#   Generate random HDF5 file of "nuclides". Parallel using simple multiprocessing
###
N_workers = int(os.getenv('SLURM_NTASKS'))
print(f"Sampling sandy with N_workers={N_workers}")

print("Beginning NJOY processing")
with Pool(N_workers) as pool:
    random_nuc = []
    for (i, nuc) in enumerate(nuclides):

        func_args = (nuc, endf_path, seed[i])
        RN = pool.apply_async(sample_nuclide_sandy, func_args)
        random_nuc.append(RN)

    for r in random_nuc:
        r.wait()

    for (i, r) in enumerate(random_nuc):
        random_nuc[i] = r.get()


###
#   Run openmc with random files
###
run_openmc(openmc_xml_dir, random_nuc,  cross_sections_xml=XS_LIB, threads = N_workers)


###
#   Store results of interest to openmc.out
###
sp = openmc.StatePoint("openmc_sim/statepoint.10.h5")

tally_TBR = sp.get_tally(name='TBR')
tally_leakage = sp.get_tally(name='leakage')

df = tally_TBR.get_pandas_dataframe()
tbr = df['mean'].sum()
tbr_std = df['std. dev.'].sum()

df_leakage = tally_leakage.get_pandas_dataframe()
leakage = df_leakage['mean'].sum()
leakage_std = df_leakage['std. dev.'].sum()

with open('openmc.out', 'w') as f:
    f.write(f"{tbr} {tbr_std} {leakage} {leakage_std}")
