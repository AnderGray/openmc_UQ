# openmc_UQ
Workflow for running nuclear data Monte Carlo UQ with OpenMC. The workflow wraps the [SANDY](https://github.com/luca-fiorito-11/sandy) and [OpenMC](https://github.com/openmc-dev/openmc) in [UncertaintyQuantification.jl](https://github.com/FriesischScott/UncertaintyQuantification.jl). The workflow generates nuclear data samples and runs openmc in the same slurm submission. Random seeds for SANDY are interpolated into an input file by `UncertaintyQuantification.jl`.

Features:
* OpenMC configuration agnostic (can use DAGMC + libmesh)
* Convenient HPC resource manager (can use multiple nodes for OpenMC simulations). Enabled through UQ.jl [HPC simulation manager](https://friesischscott.github.io/UncertaintyQuantification.jl/dev/manual/hpc)
* Sample multiple nuclides simultaneously (from available ENDF covariance files)
* Exposure to more advanced sampling algorithms from `UncertaintyQuantification.jl`
* Dimension reduction of Nuclear Data covariance matrices using SVDs

Contributions and issues are welcome.

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
7. Download nuclear data (ENDF and HDF5)
8. Install openmc_uq python package:
   - `pip install .`

## To run a basic Monte Carlo OpenMC wrap of SANDY

See `example/simple_tokamak` for example

1. Configure the `uq_inputs.json` input script

2. Run it! `julia MonteCarlo.jl uq_inputs.json`
- Or submit `julia MonteCarlo.jl uq_inputs.json` in a slurm script with a few cpus

Things to consider while configuring the input file:

- Modify `solver_exe`
- Modify your `endf_dir`
- Modify your openmc cross_sections.xml (ideally the same library as whatever ENDF library you've specified)
- Modify the `openmc_xml_dir`, should point to the directory containing your OpenMC input files
- Modify the HPC credentials
   - You can optionally change the `throttle` (max samples run simultaneous)
   - You can optionally change `ntasks` (5 cores gives a simulation time ~2 mins per simulation)
   - You can optionally change `nodes` for larger simulations
   - You can optionally set a `batchsize` (max samples submitted simultaneous) see [docs](https://friesischscott.github.io/UncertaintyQuantification.jl/dev/manual/hpc)
   - You can optionally change the number of samples. Currently set to 100 samples
- Modify the `command list`, to what needs to be loaded to run your OpenMC model
- Modify the `work_dir`, to the directory where you would like simulations to be run

This should produce the similar figures to `example/simple_tokamak/figures`. If you're not interested in computing a failure probability, the limitstate function can be ignored. You make also simulate samples without using `probability_of_failure` function. This can be done by

```julia
samples = sample(inputs)
evaluate!(ext, samples)
```
Where `inputs` and `ext` are the inputs and external model as defined in `MonteCarlo.jl`.

## To run a more advanced sampling strategy (variance reduction)

See `example/simple_tokamak_advanced` for example

1. Configure the `uq_inputs.json` input script

2. Setup the simulation `python3 setup.py uq_inputs.json`

3. Run it! `julia LatinHypercube.jl uq_inputs.json`
- Or submit `julia LatinHypercube.jl uq_inputs.json` in a slurm script with a few cpus

Other variance reduction, such as Subset Simulation are available, for evaluating design reliability, such as P(TBR < 1.05).

### What does `setup.py` do? 

It performs an SVD of all the covariances available within your specified ENDF library, stores the matrices, and computes the number of components to capture 99% of the variance of each nuclide. 

For most nuclides it is quite expensive to compute SVDs, so parallelism is also provided. In a slurm script:

```shell
SLURM_NTASKS=10 #(e.g.)
python3 setup.py uq_inputs.json
```

### Note: currently only MF33 (cross section) covariances are captured by this.
