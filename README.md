# openmc_UQ
A fairly clunky workflow for running nuclear data Monte Carlo UQ with OpenMC. The workflow wraps the [SANDY](https://github.com/luca-fiorito-11/sandy) and [OpenMC](https://github.com/openmc-dev/openmc) in [UncertaintyQuantification.jl](https://github.com/FriesischScott/UncertaintyQuantification.jl)

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

1. Modify the top path in [this line](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/run_model.py#L8), to the path to the 
2. Modify your endf path in [this line](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/run_model.py#L19)
3. Modify your openmc cross_sections.xml path in [this line](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/run_model.py#L25)
4. Change my HPC credentials to your HPC credentials (plus your desired partition etc) in [these lines](https://github.com/AnderGray/openmc_UQ/blob/4e6a457408502dbb96848ecc2bf314fc61eb2b5c/example/simple_tokamak/MonteCarlo.jl#L70)
   - You can optionally change the `throttle` (max samples run simultaneous)
   - You can optionally change `ntasks` (5 cores gives a simulation time ~2 mins per simulation)
   - You can optionally set a `batchsize` (max sample submitted simultaneous) see [docs](https://friesischscott.github.io/UncertaintyQuantification.jl/dev/manual/hpc)
 

## To run your example

In addition to the above instructions:


