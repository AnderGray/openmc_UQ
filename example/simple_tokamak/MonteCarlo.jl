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
solver_inputs            = Dict(inputs["solver_inputs"])
solver_config            = inputs["solver_config"]
reliability_qoi_name     = inputs["reliability_qoi_name"]
reliability_qoi_criteria = inputs["reliability_qoi_criteria"]
workdir                  = inputs["work_dir"]
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
println("Solver config        = $solver_config")
println("Solver inputs        = $solver_inputs")
println("Reliability QOI      = $reliability_qoi_name")
println("Reliability criteria = $reliability_qoi_criteria")
println("Work directory       = $workdir")
println("Num samples          = $num_samples")
println("Throttle             = $throttle")
println("Commands             = $command_list")
println("Slurm options        = $slurm_options")
println("Summary file         = $summary_file")
println("Results file         = $results_file")
println("******************************************************************")

############################################################################
# Specify the quantities of interest and their location in the output file

# Convenience function to retrieve data from a whitespace delimited file
function get_data(file,index)
    # Load data into an array, using whitespace as delimeter
    data = readdlm(file, ' ')
    # Retreive the data at specified index
    return data[index]
end

# Convenience function to return an extractor
# First argument to Extractor struct is itself an anonymous function
## argument base is the working directory
## function will append base with output file and extract the data
# Second argument to Extractor is the symbol name (immutable)
function openmc_extractor(index,output_file,qoi_name)
    extractor = Extractor(base -> begin
        file = joinpath(base, output_file)
	return get_data(file,index)
    end,Symbol(qoi_name))
    return extractor
end

# Loop over scores
solver_output = solver_inputs["output"]
scores = solver_inputs["scores"]
num_scores = length(scores)
extractor_list=Vector{Extractor}()
for i_score = 1:num_scores
    qoi = scores[i_score]
    # Remember, Julia doesn't use zero-indexing
    data_index = 2*i_score-1
    qoi_extractor = openmc_extractor(data_index,solver_output,qoi)
    push!(extractor_list, qoi_extractor)

    qoi_std = string(qoi,"_std")
    # Remember, Julia doesn't use zero-indexing
    std_index = 2*i_score
    qoi_extractor = openmc_extractor(std_index,solver_output,qoi_std)
    push!(extractor_list, qoi_extractor)
end
############################################################################
# Define the "solver".
solver = Solver(solver_exe, solver_config; args="")
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
source_dir = joinpath(pwd())
model = ExternalModel(
    source_dir, solver_config, extractor_list, solver; workdir=workdir, formats=numberformats, scheduler = slurm
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
nuclides=solver_inputs["nuclides"]
num_seeds=length(nuclides)
random_variable_list=Vector{RandomVariable}()
seed_string_list=Vector{String}()
for i_seed in 1:num_seeds
    seed_name=string("X","$i_seed")
    rand = RandomVariable(DiscreteUniform(1, INT64MAX), Symbol(seed_name))
    push!(random_variable_list, rand)
    seed_string=string("{{{ :","$seed_name"," }}}")
    push!(seed_string_list,seed_string)
end
# Add seeds to input dict
solver_inputs["seeds"] = seed_string_list
############################################################################
# Tell solver how much local resource to use
ntasks = parse(Int,slurm_options["ntasks"])
solver_inputs["n_cores"] = ntasks
############################################################################
# Write all inputs to a template config file
open(solver_config,"w") do f
  JSON.print(f, solver_inputs, 4)
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
println(names(samples))

println("Writing results to $results_file")
fid = h5open(results_file, "w")
fid["pf"] = pf
fid["pf_std"] = pf_std
for index = 1:num_scores
    score_name = scores[index]
    score_results = Vector{Float64}(samples[:,score_name])
    fid[score_name] = score_results
    score_std_name = string(score_name,"_std")
    score_std_results = Vector{Float64}(samples[:,score_std_name])
    fid[score_std_name] = score_std_results
end

input_slice = ["$((Symbol(:X,i)))" for i = 1:num_seeds]
fid["seeds"] = Matrix(samples[!,input_slice])
close(fid)
#############################################################################
# Generate Histograms
println("Writing histograms")
for index = 1:num_scores
    score_name = scores[index]
    histogram(samples[:,score_name], normalize=:pdf, label= "tendl2019")
    xlabel!(score_name)
    output_name=string(score_name,"_UQ.pdf")
    savefig(output_name)
    println("Wrote to $output_name")
end
println("******************************************************************")
############################################################################
