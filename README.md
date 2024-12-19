# openmc_UQ
A fairly clunky workflow for running nuclear data Monte Carlo UQ with OpenMC. The workflow wraps the [SANDY](https://github.com/luca-fiorito-11/sandy) and [OpenMC](https://github.com/openmc-dev/openmc) in [UncertaintyQuantification.jl](https://github.com/FriesischScott/UncertaintyQuantification.jl). The workflow generates nuclear data samples and runs openmc in the same slurm submission. Random seeds for SANDY are interpolated into an input file by UQ.jl.

Features:
* OpenMC configuration agnostic (can use DAGMC + libmesh)
* Convenient HPC resource manager (can use multiple nodes for OpenMC). Enabled through UQ.jl [HPC simulation manager](https://friesischscott.github.io/UncertaintyQuantification.jl/dev/manual/hpc)
* Sample multiple nuclides simultaneously (from available ENDF covariance files)

Contributions to make it less clunky welcome.

## Installation

1. Install [OpenMC](https://github.com/openmc-dev/openmc)
2. Install [NJOY](https://github.com/njoy/NJOY21)
3. Install [SANDY](https://github.com/luca-fiorito-11/sandy)
   - `pip install sandy`
   - Or from [source](https://github.com/luca-fiorito-11/sandy/blob/develop/INSTALL.md#installing-sandy-from-source)
5. Install [Julia](https://julialang.org/downloads/). Recommended to use the Juliaup installation manager
6. Install [UncertaintyQuantification.jl](https://github.com/FriesischScott/UncertaintyQuantification.jl) and other julia packages
   - `julia -e 'using Pkg; Pkg.add("UncertaintyQuantification")'`
   - `julia -e 'using Pkg; Pkg.add("DelimitedFiles")'`
   - `julia -e 'using Pkg; Pkg.add("HDF5")'`
   - `julia -e 'using Pkg; Pkg.add("StatsBase")'`
   - `julia -e 'using Pkg; Pkg.add("Plots")'`
8. Download nuclear data (ENDF and HDF5)
 
## To run my example

1. Modify your top path in [this line](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/run_model.py#L8) (the top directory of this repo)
2. Modify your endf path in [this line](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/run_model.py#L19)
3. Modify your openmc cross_sections.xml path in [this line](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/run_model.py#L25)
4. Change my HPC credentials to your HPC credentials (plus your desired partition etc) in [these lines](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/MonteCarlo.jl#L70)
   - You can optionally change the `throttle` (max samples run simultaneous)
   - You can optionally change `ntasks` (5 cores gives a simulation time ~2 mins per simulation)
   - You can optionally set a `batchsize` (max samples submitted simultaneous) see [docs](https://friesischscott.github.io/UncertaintyQuantification.jl/dev/manual/hpc)
   - You can optionally change the number of samples [here](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/MonteCarlo.jl#L96). Currently set to 1000 samples
     
5. Modify the "extras" that will be injected into slurm script to run your model, in [this line](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/MonteCarlo.jl#L82)
    - For me, this was loading openmc, activating my python environment, and loading NJOY
6. Run it! From inside `examples/simple_tokamak`
    - `julia`
    - `include("MonteCarlo.jl")`
    - Or submit `julia -e 'include("MonteCarlo.jl")'` in a slurm script with 1 task
    
This should produce the similar figures to `example/simple_tokamak/figures`. If you're not interested in computing a failure probability, the limitstate function can be ignored. You make also simulate samples without using `probability_of_failure` function. This can be done by

```julia
samples = sample(inputs)
evaluate!(ext, samples)
```
Where `inputs` and `ext` are the inputs and external model as defined in `MonteCarlo.jl`. The tallies of interest must then be extracted manually however.

## To run your example

In addition to the above instructions:

1. Make a new folder in `examples`, and copy `MonteCarlo.jl` and `run_model.py`
2. Change [this directory](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/run_model.py#L22) to the directory of your model (.xml files, .h5m files, ...)
3. Change the list of nuclides you wish to sample [here](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/run_model.py#L17)
4. Output your tally of interest [in lines 65-79](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/run_model.py#L65)
5. Change the extractor to your tally of interest [here](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/MonteCarlo.jl#L33)
   - Currently only scalar tallies can be "extracted", but you can generate samples of any tally (e.g. spectra, meshes), and extract them yourself.
     
6. Change the `nodes` and `ntasks` to your desired resources, [here](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/MonteCarlo.jl#L74)
   - If you increase the `nodes` to more than 1, you must match it [here](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/src/run_openmc.py#L64) 😅
       
8. 😅 Change the length of the 'seeds' [here](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/run_model.py#L31) to the length of the nuclides
   - If `nuclides = ["Fe54", "Fe56"]`, then `seeds = [{{{   :X1   }}}, {{{   :X2   }}}]`
   - If `len(nuclides) == 50`, then `seeds = [{{{   :X1   }}}, {{{   :X2   }}}, ... , {{{   :X50   }}}]`
9. 😅 Do the same for the RandomVariables [here](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/MonteCarlo.jl#L6) and [here](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/MonteCarlo.jl#L17)
    
10. Run it!

The "😅" are locations where the workflow can be clearly improved (among other places)
