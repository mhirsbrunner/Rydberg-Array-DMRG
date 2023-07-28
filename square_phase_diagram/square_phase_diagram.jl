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

@passobj 1 workers() n_threads
@passobj 1 workers() config_dir
@passobj 1 workers() data_dir
@passobj 1 workers() outputlevel

@everywhere BLAS.set_num_threads(n_threads)

if parsed_args["write_dir"] == nothing
    @everywhere function func(a)
        RydbergDMRG.run_phase_diagram_point(a[1], a[2], config_dir, data_dir; outputlevel=outputlevel)
    end
else  
    write_dir = parsed_args["write_dir"]
    @passobj 1 workers() write_dir

    @everywhere function func(a)
        RydbergDMRG.run_phase_diagram_point(a[1], a[2], config_dir, data_dir; outputlevel=outputlevel, write_dir=write_dir)
    end
end

phase_space_config = JSON3.read(read(joinpath(config_dir, "phase_space_config.json")))

rb_over_a_ax = phase_space_config["rb_over_a_min"]:phase_space_config["rb_over_a_step"]:phase_space_config["rb_over_a_max"]
delta_over_omega_ax = phase_space_config["delta_over_omega_min"]:phase_space_config["delta_over_omega_step"]:phase_space_config["delta_over_omega_max"]

Parallelism.robust_pmap(func, Iterators.product(rb_over_a_ax, delta_over_omega_ax))