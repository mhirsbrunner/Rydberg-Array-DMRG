using ArgParse
using Distributed
using ThreadPinning

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
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

n_threads = parsed_args["n_threads"]
config_dir = parsed_args["config_dir"]
data_dir = parsed_args["data_dir"]
outputlevel = parsed_args["outputlevel"]

phase_space_config = JSON3.read(read(joinpath(config_dir, "local_detuning_phase_space_config.json")))

rb = phase_space_config["rb"]
delta_ax = phase_space_config["delta_min"]:phase_space_config["delta_step"]:phase_space_config["delta_max"]
delta_local_ax = phase_space_config["delta_local_min"]:phase_space_config["delta_local_step"]:phase_space_config["delta_local_max"]

@passobj 1 workers() n_threads
@passobj 1 workers() config_dir
@passobj 1 workers() data_dir
@passobj 1 workers() outputlevel
@passobj 1 workers() rb

@everywhere BLAS.set_num_threads(n_threads)

if parsed_args["write_dir"] == nothing
    @everywhere function func(a)
        RydbergDMRG.run_local_detuning_phase_diagram_point(rb, a[1], a[2], config_dir, data_dir; outputlevel=outputlevel)
    end
else  
    write_dir = parsed_args["write_dir"]
    @passobj 1 workers() write_dir

    @everywhere function func(a)
        RydbergDMRG.run_local_detuning_phase_diagram_point(rb, a[1], a[2], config_dir, data_dir; outputlevel=outputlevel, write_dir=write_dir)
    end
end

Parallelism.robust_pmap(func, Iterators.product(delta_ax, delta_local_ax))