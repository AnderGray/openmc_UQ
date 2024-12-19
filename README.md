# openmc_UQ
A fairly clunky workflow for running nuclear data Monte Carlo UQ with OpenMC. The workflow wraps the [SANDY](https://github.com/luca-fiorito-11/sandy) and [OpenMC](https://github.com/openmc-dev/openmc) in [UncertaintyQuantification.jl](https://github.com/FriesischScott/UncertaintyQuantification.jl)

Features:
* OpenMC configuration agnostic (can use DAGMC + libmesh)
* Sample multiple nuclides simultaneously (from available ENDF covariance files)
* Convenient HPC resource manager (can use multiple nodes for OpenMC). Enabled through UQ.jl [HPC simulation manager](https://friesischscott.github.io/UncertaintyQuantification.jl/dev/manual/hpc)


## Installation

1. Install [OpenMC](https://github.com/openmc-dev/openmc)
2. Install [NJOY](https://github.com/njoy/NJOY21)
3. Install [SANDY](https://github.com/luca-fiorito-11/sandy)
   - `pip install sandy`
   - Or from [source](https://github.com/luca-fiorito-11/sandy/blob/develop/INSTALL.md#installing-sandy-from-source)
5. Install [Julia](https://julialang.org/downloads/). Recommended to use the Juliaup installation manager
6. Install [UncertaintyQuantification.jl](https://github.com/FriesischScott/UncertaintyQuantification.jl) and other julia packages
   - `julia -e 'using Pkg; Pkg.add("UncertaintyQuantification")'`
   - `julia -e 'using Pkg; Pkg.add("HDF5")'`
   - `julia -e 'using Pkg; Pkg.add("StatsBase")'`
   - `julia -e 'using Pkg; Pkg.add("Plots")'`
8. Download nuclear data (ENDF and HDF5)
 
## To run my example

## To run your example
