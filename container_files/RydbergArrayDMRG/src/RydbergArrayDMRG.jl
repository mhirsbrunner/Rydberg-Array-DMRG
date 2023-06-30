module RydbergArrayDMRG

using JSON3
using ITensors

C6 = 2π * 862690; # MHz μm^6

function rydberg_ham_1d(n_sites, omega, delta, a, n_nns=4)
    if n_nns >= n_sites
        n_nns = n_sites - 1
    end

    sites = siteinds("S=1/2", n_sites)

    os = OpSum()

    for j=1:n_sites
        os += omega, "Sx", j
        os += -delta, "ProjUp", j
    end

    for k=1:n_nns
        for j=1:n_sites - k
            r = a * k
            V = C6 / (r ^ 6)
            os += V, "ProjUp", j, "ProjUp", j + k
        end
    end

    hamiltonian = MPO(os, sites)

    return hamiltonian, sites
end

function rydberg_ham_2d(n_x, n_y, omega, delta, a, yperiodic=true, n_nns=4)
    if n_nns >= min(n_x, n_y)
        n_nns = min(n_x, n_y) - 1
    end
    
    n_sites = n_x * n_y

    sites = siteinds("S=1/2", n_sites)

    lattice = square_lattice(n_x, n_y; yperiodic=yperiodic)
    os = OpSum()

    for j=1:n_sites
        os += omega, "Sx", j
        os += -delta, "ProjUp", j
    end

    for b=1:n_nns
        for j=1:n_sites - b
            r = a * b
            V = C6 / (r ^ 6)
            os += V, "ProjUp", j, "ProjUp", j + b
        end
    end

    hamiltonian = MPO(os, sites)

    return hamiltonian, sites
end

function rydberg_density(psi)
    return expect(psi, "ProjUp")
end

function dmrg_run_1d_chain(rb_over_a, delta_over_omega, config_dir="config", data_dir="data_dir"; outputlevel=0, write_when_maxdim_exceeds=5000, write_dir="SCRATCH")
    # Hamiltonian specification
    ham_config = JSON3.read(read(joinpath(config_dir, "1d_ham_config.json")))
    n_sites = ham_config["n_sites"]
    omega = 2π * ham_config["omega"]
    delta = delta_over_omega * omega

    Rb = (C6 / omega) ^ (1 / 6)
    a = Rb / rb_over_a

    hamiltonian_params = Dict("n_sites"=>n_sites, "omega"=>omega, "delta"=>delta, "a"=>a)
    
    if haskey(ham_config, "n_nns")
        n_nns = hamiltonian_params["n_nns"]
        ham, sites = rydberg_ham_1d(n_sites, omega, delta, a, n_nns)
        hamiltonian_params["n_nns"] = n_nns
    else
        ham, sites = rydberg_ham_1d(n_sites, omega, delta, a)
    end

    # DMRG specifications
    dmrg_params = JSON3.read(read(joinpath(config_dir, "1d_dmrg_config.json")))
    nsweeps = dmrg_params["nsweeps"]
    maxdim = dmrg_params["maxdim"]
    cutoff = dmrg_params["cutoff"]
    krylovdim = dmrg_params["krylovdim"]

    # Run DMRG
    psi0 = randomMPS(sites, maxdim[1])
    energy_observer = DMRGObserver(energy_tol=dmrg_params["energy_tol"])

    write_path = joinpath(write_dir, "$(rb_over_a)_$(delta_over_omega)")
    mkpath(write_path)
    energy, psi = dmrg(ham, psi0; nsweeps, maxdim, cutoff, eigsolve_krylovdim=krylovdim, observer=energy_observer, outputlevel=outputlevel, write_when_maxdim_exceeds=write_when_maxdim_exceeds, write_path=write_path)
    rm(write_path, recursive=true)

    # Save results
    params = Dict("hamiltonian_params" => hamiltonian_params, "dmrg_params" => dmrg_params)
    result = Dict("energy" => energy, "rydberg_density" => rydberg_density(psi))

    mkpath(data_dir)
    JSON3.write(joinpath(data_dir, "$(rb_over_a)-$(delta_over_omega)" * ".json"), Dict("result" => result, "params" => params))
end

end # module RydbergArrayDMRG
