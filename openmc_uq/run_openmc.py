import os
from pathlib import Path
from shutil import copyfile
import openmc.data
from .utils import replace_nuclide_tally, replace_nuclide_material

def run_openmc(openmc_xml_dir, random_nuclides, cross_sections_xml,
               threads = 1,
               run_dir="openmc_sim"):

    # ==============================================================================
    # Make openmc sim directory

    working_dir = os.getcwd()

    print()
    print(" *************** Making sim folder ***************")
    print()

    out_dir = Path(run_dir).resolve()
    out_dir.mkdir(exist_ok=True)

    openmc_xml_dir = Path(openmc_xml_dir).resolve()
    input_files  = ["materials.xml", "settings.xml", "tallies.xml", "geometry.xml"]

    for file in input_files:
        copyfile(openmc_xml_dir/file, out_dir/file)    

    for file in os.listdir(openmc_xml_dir):
        if file.endswith(".h5m"):
            copyfile(openmc_xml_dir/file, out_dir/file)
        if file.endswith(".e"):
            copyfile(openmc_xml_dir/file, out_dir/file)    

    os.chdir(out_dir)

    # ==============================================================================
    # Create xml library

    print()
    print(" *************** Creating random library ***************")
    print()

    lib = openmc.data.DataLibrary()
    lib = lib.from_xml(cross_sections_xml)  # Gets current

    for nuc in random_nuclides:
        lib.register_file(nuc.path)

    post = out_dir / "cross_sections_rand.xml"
    lib.export_to_xml(post)

    # ==============================================================================
    # Change openmc inputs

    [replace_nuclide_material(nuc.nuclide, nuc.perturbed_name) for nuc in random_nuclides]
    [replace_nuclide_tally(nuc.nuclide, nuc.perturbed_name) for nuc in random_nuclides]

    materials   = openmc.Materials.from_xml("materials.xml")
    materials.cross_sections = str(post)
    materials.export_to_xml()

    #output = openmc.run(threads = threads)
    openmc_command = f"mpirun -np 1 -ppn 1 --bind-to none openmc -s {threads}"
    os.system(openmc_command)

    os.chdir(working_dir)

    #return output


