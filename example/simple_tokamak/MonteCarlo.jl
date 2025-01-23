using UncertaintyQuantification, DelimitedFiles, HDF5, StatsBase, Plots

############################################################################
## Things to un-hardcode
## Number of random variables for nuclides
## Workdir
## Source directory - make it possible to run from elsewhere
## source file
## HPC credentials and settings
## For loop over QOI
## Openmc out file
## Number of samples
## throttle (max samples run simultaneous)
## You can optionally set a batchsize (max samples submitted simultaneous)
## see https://friesischscott.github.io/UncertaintyQuantification.jl/dev/manual/hpc

## Histograms, output files

## Quality of life improvments
# Character limit
# Refactor into functions

# These files will be rendered through Mustach.jl and have values injected
sourcedir = joinpath(pwd())
sourcefile = "run_model.py"

# UQ.jl will create subfolders in here to run the solver and store the results
workdir = joinpath(pwd(), "run_dir")

# HPC options
Slurm_options = Dict(
"account"=>"UKAEA-AP001-CPU",
"partition"=>"icelake",
"job-name"=>"openmc_UQ",
"nodes"=>"1",
"ntasks" =>"5",
"time"=>"00:20:00"
)

# Extra commands
command_list=["source ~/work/opt/openmc_intel_helen/profiles/openmc_profile", "source ~/.virtualenvs/openmc/bin/activate", "module load njoy/21", "export NJOY=njoy"]
############################################################################

###
#   Random seeds for each nuclide. Must be the same number as nuclides. The big number is just typemax(Int64)
###
X1 = RandomVariable(DiscreteUniform(1, 9223372036854775807), :X1)
X2 = RandomVariable(DiscreteUniform(1, 9223372036854775807), :X2)
X3 = RandomVariable(DiscreteUniform(1, 9223372036854775807), :X3)
X4 = RandomVariable(DiscreteUniform(1, 9223372036854775807), :X4)
X5 = RandomVariable(DiscreteUniform(1, 9223372036854775807), :X5)
X6 = RandomVariable(DiscreteUniform(1, 9223372036854775807), :X6)
X7 = RandomVariable(DiscreteUniform(1, 9223372036854775807), :X7)
X8 = RandomVariable(DiscreteUniform(1, 9223372036854775807), :X8)
X9 = RandomVariable(DiscreteUniform(1, 9223372036854775807), :X9)
X10 = RandomVariable(DiscreteUniform(1, 9223372036854775807), :X10)

inputs = [X1, X2, X3, X4, X5, X6, X7, X8, X9, X10]

# Dictionary to map format Strings (Formatting.jl) to variables
numberformats = Dict(:* => "d")


# Read output files. Specify the quantities of interest and their location in the output file
# An extractor is based the working directory for the current sample
TBR = Extractor(base -> begin
    file = joinpath(base, "openmc.out")
    data = readdlm(file, ' ')

    return data[1]
end, :TBR)

TBR_std = Extractor(base -> begin
    file = joinpath(base, "openmc.out")
    data = readdlm(file, ' ')

    return data[2]
end, :TBR_std)

leakage = Extractor(base -> begin
    file = joinpath(base, "openmc.out")
    data = readdlm(file, ' ')

    return data[3]
end, :leakage)

leakage_std = Extractor(base -> begin
    file = joinpath(base, "openmc.out")
    data = readdlm(file, ' ')

    return data[4]
end, :leakage_std)

# Your "solver". Can be anything called from the command line
openmc = Solver(
    "python3", # path to OpenSees binary
    sourcefile;
    args="", # (optional) extra arguments passed to the solver
)

slurm = SlurmInterface(
    Slurm_options;
    throttle=200,
    extras=command_list,
    )

ext = ExternalModel(
    sourcedir, sourcefile, [leakage, leakage_std, TBR, TBR_std], openmc; workdir=workdir, formats=numberformats, scheduler = slurm
)


## Limistate function for reliability analysis. Here computing Prob(TBR <= 1.05)
function limitstate(df)
    return  reduce(vcat, df.TBR) .- 1.05
end


sim = MonteCarlo(1000)      # Change number of samples

println("Starting Monte Carlo simulation")
@time pf, pf_std, samples = probability_of_failure(ext, limitstate, inputs, sim)       # Run UQ

println("Probability of failure MC: $pf")
println("pf_std: $pf_std")

file = open("probability_of_failure_MC.txt", "w")
write(file, "pf = $pf \n pf_std = $pf_std")
close(file)

# Save results to HDF5
input_slice = ["$((Symbol(:X,i)))" for i = 1:length(inputs)]

fid = h5open("result_mc", "w")
fid["pf"] = pf
fid["pf_std"] = pf_std
fid["TBR"] = samples.TBR
fid["leakage"] = samples.leakage
fid["inputs"] = Matrix(samples[!,input_slice])
close(fid)

# Print stats
TBR = samples.TBR #hide
TBR_mean = mean(TBR) #hide
TBR_std  = std(TBR) #hide
lower_quantile = quantile(TBR, 0.025) #hide
upper_quantile = quantile(TBR, 0.975) #hide
println("TBR mean: $TBR_mean, TBR std: $TBR_std, TBR 95%: [$lower_quantile, $upper_quantile]") #hide

# Generate Histograms
histogram(samples.TBR, normalize=:pdf, label= "tendl2019")
xlabel!("TBR")
savefig("TBR_UQ.pdf")

histogram(samples.leakage, normalize=:pdf, label= "tendl2019")
xlabel!("leakage")
savefig("leakage_UQ.pdf")
