using UncertaintyQuantification, DelimitedFiles, HDF5, StatsBase, Plots, JSON

############################################################################
## Things to un-hardcode
## Number of random variables for nuclides
## # You can optionally set a batchsize (max samples submitted simultaneous)
## see https://friesischscott.github.io/UncertaintyQuantification.jl/dev/manual/hpc
############################################################################
## Quality of life improvments
# Character limit
# Refactor into functions
############################################################################
# Process command line arguments
input_file=""
try
    global input_file = ARGS[1]
catch e
    error("Please provide an input JSON file. Usage:\n","julia ",PROGRAM_FILE, " input.json")
end
#############################################################################
# Get inputs from json file
println("Reading inputs from: ", input_file)
inputs = JSON.parsefile(input_file)

# Set local variables
solver_exe               = inputs["solver_exe"]
solver_input             = inputs["solver_input"]
solver_args              = inputs["solver_args"]
solver_output            = inputs["solver_output"]
solver_dir               = inputs["solver_dir"]
solver_qoi_names         = Vector{String}(inputs["solver_qoi_names"])
reliability_qoi_name     = inputs["reliability_qoi_name"]
reliability_qoi_criteria = inputs["reliability_qoi_criteria"]
workdir                  = inputs["work_dir"]
num_seeds                = Int(inputs["num_seeds"])
num_samples              = Int(inputs["num_samples"])
throttle                 = Int(inputs["throttle"])
command_list             = Vector{String}(inputs["command_list"])
slurm_options            = Dict{String,String}(inputs["slurm_options"])
summary_file             = inputs["summary_file"]
results_file             = inputs["results_file"]

# Logging
println("******************************************************************")
println("Run configuration:")
println("******************************************************************")
println("Solver executable    = $solver_exe")
println("Solver input         = $solver_input")
println("Solver args          = $solver_args")
println("Solver output        = $solver_output")
println("Solver directory     = $solver_dir")
println("Solver QOI names     = $solver_qoi_names")
println("Reliability QOI      = $reliability_qoi_name")
println("Reliability criteria = $reliability_qoi_criteria")
println("Work directory       = $workdir")
println("Num seeds            = $num_seeds")
println("Num samples          = $num_samples")
println("Throttle             = $throttle")
println("Commands             = $command_list")
println("Slurm options        = $slurm_options")
println("Summary file         = $summary_file")
println("Results file         = $results_file")
println("******************************************************************")
############################################################################
# Specify the quantities of interest and their location in the output file

# Function to return an extractor
# First argument to Extractor struct is itself an anonymous function
## argument base is the working directory
## function will append base with output file and extract the data
# Second argument to Extractor is the symbol name (immutable once changed
function openmc_extractor(index,output_file, qoi_name)
    extractor = Extractor(base -> begin
        file = joinpath(base, output_file)
        data = readdlm(file, ' ')
	return data[index]
    end,Symbol(qoi_name))
    return extractor
end

# NB Julia array indexing starts at 1
num_qoi = length(solver_qoi_names)
extractor_list=Vector{Extractor}()
for index = 1:num_qoi
    qoi = solver_qoi_names[index]
    qoi_extractor = openmc_extractor(index,solver_output,qoi)
    push!(extractor_list, qoi_extractor)
end
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
function limitstate(df)
    return reduce(vcat, df[:,reliability_qoi_name]).-reliability_qoi_criteria
end
############################################################################
# Define sampling method (Monte Carlo)
sampling = MonteCarlo(num_samples)
############################################################################
# Define random variables corresponding to seeds for each nuclide.
INT64MAX=typemax(Int64)
random_variable_list=Vector{RandomVariable}()
for i_seed in 1:num_seeds
    seed_name=string("X","$i_seed")
    rand = RandomVariable(DiscreteUniform(1, INT64MAX), Symbol(seed_name))
    push!(random_variable_list, rand)
end
############################################################################
# Peform the Monte Carlo UQ
println("Starting Monte Carlo simulation")
@time pf, pf_std, samples = probability_of_failure(model, limitstate, random_variable_list, sampling)
println("Monte Carlo simulation complete")
############################################################################
# Extract reliability results
qoi_results = samples[:,reliability_qoi_name]
qoi_mean = mean(qoi_results)
qoi_std  = std(qoi_results)
lower_quantile = quantile(qoi_results, 0.025)
upper_quantile = quantile(qoi_results, 0.975)
############################################################################
# Logging Summary data
println("******************************************************************")
println("Summary:")
println("******************************************************************")
println("Probability of failure: $pf")
println("Probability of failure standard deviation: $pf_std")
println("QOI mean: $qoi_mean")
println("QOI std: $qoi_std")
println("QOI 95% confidence interval: [$lower_quantile, $upper_quantile]")
println("******************************************************************")
############################################################################
# Write summary to text file
println("Writing summary to $summary_file")
file = open(summary_file, "w")
write(file, "pf = $pf \n pf_std = $pf_std")
close(file)
############################################################################
# Save results to HDF5
println("Writing results to $results_file")
input_slice = ["$((Symbol(:X,i)))" for i = 1:length(random_variable_list)]
fid = h5open(results_file, "w")
fid["pf"] = pf
fid["pf_std"] = pf_std
for index = 1:num_qoi
    qoi = solver_qoi_names[index]
    fid[qoi] = samples[:,qoi]
end
fid["input_seeds"] = Matrix(samples[!,input_slice])
close(fid)
#############################################################################
# Generate Histograms
println("Writing histograms")
for index = 1:num_qoi
    qoi = solver_qoi_names[index]
    histogram(samples[:,qoi], normalize=:pdf, label= "tendl2019")
    xlabel!(qoi)
    output_name=string(qoi,"_UQ.pdf")
    savefig(output_name)
    println("Wrote to $output_name")
end
println("******************************************************************")
############################################################################
