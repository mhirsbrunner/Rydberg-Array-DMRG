using ArgParse
using Distributed
using ThreadPinning

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "n_y"
            arg_type = Int
        "n_workers"
            arg_type = Int
        "n_threads"
            arg_type = Int
        "config_dir"
            arg_type = String
        "data_dir"
            arg_type = String
        "--write_dir"
            arg_type = String
        "--outputlevel"
            arg_type = Int
            default = 0
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

addprocs(parsed_args["n_workers"]; enable_threaded_blas=true)
pinthreads(:compact)

@everywhere begin
    using ParallelDataTransfer
    using LinearAlgebra
    ITensors.Strided.disable_threads()
    include("/global/homes/m/mhirsbru/Rydberg-Array-DMRG/RydbergDMRG/RydbergDMRG.jl")
end

@everywhere using .RydbergDMRG

n_y = parsed_args["n_y"]
n_threads = parsed_args["n_threads"]
config_dir = parsed_args["config_dir"]
data_dir = parsed_args["data_dir"]
outputlevel = parsed_args["outputlevel"]

@passobj 1 workers() n_y
@passobj 1 workers() n_threads
@passobj 1 workers() config_dir
@passobj 1 workers() data_dir
@passobj 1 workers() outputlevel

@everywhere BLAS.set_num_threads(n_threads)

phase_space_config = JSON3.read(read(joinpath(config_dir, "finite_size_config.json")))

rb = phase_space_config["rb"]
op_name = phase_space_config["op_name"]
@passobj 1 workers() rb
@passobj 1 workers() op_name

delta_ax = phase_space_config["delta_min"]:phase_space_config["delta_step"]:phase_space_config["delta_max"]

if parsed_args["write_dir"] == nothing
    @everywhere function func(input)
        RydbergDMRG.calculate_finite_size_scaling_data(n_y, rb, input, config_dir, data_dir; outputlevel=outputlevel, op_name=op_name)
    end
else  
    write_dir = parsed_args["write_dir"]
    @passobj 1 workers() write_dir

    @everywhere function func(input)
        RydbergDMRG.calculate_finite_size_scaling_data(n_y, rb, input, config_dir, data_dir; outputlevel=outputlevel, write_dir=write_dir, op_name=op_name)
    end
end

Parallelism.robust_pmap(func, delta_ax)