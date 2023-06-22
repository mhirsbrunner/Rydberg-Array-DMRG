using JSON3
using ITensors
using ArgParse
using LinearAlgebra

BLAS.set_num_threads(Threads.nthreads())
ITensors.Strided.disable_threads()

include("rydberg_arrays.jl")
C6 = 2π * 862690; # MHz μm^6

function parse_commandline()
    s = ArgParseSettings()
    
    @add_arg_table s begin

        "rb_over_a"
            arg_type = Float64
            required = true

        "delta_over_omega"
            arg_type = Float64
            required = true

        "--data_dir"
            arg_type = String
            default = "phase_diagram_data"
        
        "--fname"
            arg_type = String
            default = "delta_sweep_data"

        "--output_level"
            arg_type = Int
            default = 1
    end

    return parse_args(s)

end

function main()
    parsed_args = parse_commandline()

    # Hamiltonian specification
    ham_config = JSON3.read(read("1d_ham_config.json"))
    n_sites = ham_config["n_sites"]
    omega = 2π * ham_config["omega"]
    delta = parsed_args["delta_over_omega"] * omega

    Rb = (C6 / omega) ^ (1 / 6)
    a = Rb / parsed_args["rb_over_a"]

    hamiltonian_params = Dict("n_sites"=>n_sites, "omega"=>omega, "delta"=>delta, "a"=>a)
    
    if haskey(ham_config, "n_nns")
        n_nns = hamiltonian_params["n_nns"]
        ham, sites = rydberg_ham_1d(n_sites, omega, delta, a, n_nns)
        hamiltonian_params["n_nns"] = n_nns
    else
        ham, sites = rydberg_ham_1d(n_sites, omega, delta, a)
    end

    # DMRG specifications
    dmrg_params = JSON3.read(read("1d_dmrg_config.json"))
    nsweeps = dmrg_params["nsweeps"]
    maxdim = dmrg_params["maxdim"]
    cutoff = dmrg_params["cutoff"]
    krylovdim = dmrg_params["krylovdim"]

    # Run DMRG
    psi0 = randomMPS(sites, maxdim)
    energy_observer = DMRGObserver(energy_tol=dmrg_params["energy_tol"])

    energy, psi = dmrg(ham, psi0; nsweeps, maxdim, cutoff, eigsolve_krylovdim=krylovdim, observer=energy_observer, outputlevel=parsed_args["output_level"])

    # Save results
    params = Dict("hamiltonian_params" => hamiltonian_params, "dmrg_params" => dmrg_params)
    result = Dict("energy" => energy, "rydberg_density" => rydberg_density(psi))

    JSON3.write(joinpath(parsed_args["data_dir"], parsed_args["fname"] * ".json"), Dict("result" => result, "params" => params))
end

main()
