# %% Imports
import os
import openmc.data
from pathlib import Path
from .utils import get_nuclide_paths

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
