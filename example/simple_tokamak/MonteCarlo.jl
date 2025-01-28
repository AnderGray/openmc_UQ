using UncertaintyQuantification, DelimitedFiles, HDF5, StatsBase, Plots, JSON

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
## # You can optionally set a batchsize (max samples submitted simultaneous)
## see https://friesischscott.github.io/UncertaintyQuantification.jl/dev/manual/hpc
## Histograms, output files
############################################################################
## Quality of life improvments
# Character limit
# Refactor into functions
############################################################################
# Get inputs from json file
# TODO make this an argument
filename="uq_inputs.json"
inputs = JSON.parsefile(filename)

# Set local variables
solver_exe    = inputs["solver_exe"]
solver_input  = inputs["solver_input"]
solver_args   = inputs["solver_args"]
solver_output = inputs["solver_output"]
solver_dir    = inputs["solver_dir"]
workdir       = inputs["work_dir"]
num_samples   = Int(inputs["num_samples"])
throttle      = Int(inputs["throttle"])
command_list  = Vector{String}(inputs["command_list"])
slurm_options = Dict{String,String}(inputs["slurm_options"])

# Logging
println("************************************************")
println("Run configuration:")
println("************************************************")
println("Solver executable = $solver_exe")
println("Solver input      = $solver_input")
println("Solver args       = $solver_args")
println("Solver output     = $solver_output")
println("Solver directory  = $solver_dir")
println("Work directory    = $workdir")
println("Num samples       = $num_samples")
println("Throttle          = $throttle")
println("Commands          = $command_list")
println("Slurm options     = $slurm_options")
println("************************************************")

############################################################################
# Specify the quantities of interest and their location in the output file
# An extractor is based the working directory for the current sample
# TODO: make into a loop
TBR = Extractor(base -> begin
    file = joinpath(base, solver_output)
    data = readdlm(file, ' ')

    return data[1]
end, :TBR)

TBR_std = Extractor(base -> begin
    file = joinpath(base, solver_output)
    data = readdlm(file, ' ')

    return data[2]
end, :TBR_std)

leakage = Extractor(base -> begin
    file = joinpath(base, solver_output)
    data = readdlm(file, ' ')

    return data[3]
end, :leakage)

leakage_std = Extractor(base -> begin
    file = joinpath(base, solver_output)
    data = readdlm(file, ' ')

    return data[4]
end, :leakage_std)

extractor_list =  [leakage, leakage_std, TBR, TBR_std]

############################################################################
# Define the "solver".
# solver_exe can be anything called from the command line
# solver_args: (optional) extra arguments passed to the solver
solver = Solver(solver_exe, solver_input; args=solver_args)
############################################################################
# Define the job scheduler
slurm = SlurmInterface(slurm_options; throttle=throttle,extras=command_list)
############################################################################
# Define the number format
# Dictionary to map format Strings (Formatting.jl) to variables
numberformats = Dict(:* => "d")
############################################################################
# Put everything together and define (external) model
# workdir: create subfolders in here to run the solver and store the results
model = ExternalModel(
    solver_dir, solver_input, extractor_list, solver; workdir=workdir, formats=numberformats, scheduler = slurm
)
############################################################################
# Define Limistate function for reliability analysis.
# Here, computing Prob(TBR <= 1.05)
function limitstate(df)
    return  reduce(vcat, df.TBR) .- 1.05
end
############################################################################
# Define sampling method (Monte Carlo)
sampling = MonteCarlo(num_samples)
############################################################################
#
# Random seeds for each nuclide. Must be the same number as nuclides.
# The big number is just typemax(Int64)
#
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

input_seeds = [X1, X2, X3, X4, X5, X6, X7, X8, X9, X10]
############################################################################
# Peform the Monte Carlo UQ
println("Starting Monte Carlo simulation")
@time pf, pf_std, samples = probability_of_failure(model, limitstate, input_seeds, sampling)
println("Monte Carlo simulation complete")

############################################################################
# Logging Summary data
println("Probability of failure MC: $pf")
println("pf_std: $pf_std")
# Print stats
TBR = samples.TBR #hide
TBR_mean = mean(TBR) #hide
TBR_std  = std(TBR) #hide
lower_quantile = quantile(TBR, 0.025) #hide
upper_quantile = quantile(TBR, 0.975) #hide
println("TBR mean: $TBR_mean, TBR std: $TBR_std, TBR 95%: [$lower_quantile, $upper_quantile]") #hide
############################################################################
# Write summary to text file
println("Writing summary")
summary_file_name = "probability_of_failure_MC.txt"
file = open(summary_file_name, "w")
write(file, "pf = $pf \n pf_std = $pf_std")
close(file)

############################################################################
# Save results to HDF5
println("Writing results")
results_file = "result_mc"
input_slice = ["$((Symbol(:X,i)))" for i = 1:length(input_seeds)]
fid = h5open(results_file, "w")
fid["pf"] = pf
fid["pf_std"] = pf_std
fid["TBR"] = samples.TBR
fid["leakage"] = samples.leakage
fid["input_seeds"] = Matrix(samples[!,input_slice])
close(fid)

#############################################################################
# Generate Histograms
println("Writing Histograms")
histogram(samples.TBR, normalize=:pdf, label= "tendl2019")
xlabel!("TBR")
savefig("TBR_UQ.pdf")

histogram(samples.leakage, normalize=:pdf, label= "tendl2019")
xlabel!("leakage")
savefig("leakage_UQ.pdf")
println("Writing Histograms")
############################################################################