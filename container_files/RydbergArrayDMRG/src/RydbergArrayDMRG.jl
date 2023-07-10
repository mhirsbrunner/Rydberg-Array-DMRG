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

function generate_sites_square_lattice(n_x, n_y)
    sites = []

    y_ind = 1
    for ii in 1:n_x
        for jj in 1:n_y
            push!(sites, [ii, y_ind])

            if isodd(ii)
                y_ind += 1
            else
                y_ind -= 1
            end
        end
        if isodd(ii)
            y_ind -= 1
        else
            y_ind += 1
        end
    end

    return sites
end

function site_to_index_square_lattice(site, n_x, n_y)
    if isodd(site[1])
        return n_y * (site[1] - 1) + site[2]
    else
        return n_y * (site[1] - 1) + 5 - site[2]
    end
end

function rydberg_ham_square_lattice(n_x, n_y, omega, delta, a; yperiodic=true, rr_cutoff=0.05)    
    n_sites = n_x * n_y
    
    lattice_sites = generate_sites_square_lattice(n_x, n_y)

    os = OpSum()

    for ii = 1:n_sites
        os += omega, "Sx", ii
        os += -delta, "ProjUp", ii
    end

    for ind_1 = 1:n_sites
        for ind_2 = ind_1 + 1:n_sites
            site_1 = lattice_sites[ind_1]
            site_2 = lattice_sites[ind_2]

            x_disp = site_1[1] - site_2[1]
            temp_y_disp = site_1[2] - site_2[2]

            if yperiodic
                y_disp = min(abs(temp_y_disp), n_y - abs(temp_y_disp))
            else
                y_disp = temp_y_disp
            end

            r = a * sqrt(x_disp ^ 2 + y_disp ^ 2)
            V = C6 / (r ^ 6)

            if V >= rr_cutoff
                os += V, "ProjUp", ind_1, "ProjUp", ind_2
            end
        end
    end

    sites = siteinds("S=1/2", n_sites)

    hamiltonian = MPO(os, sites)

    return hamiltonian, sites
end

function rydberg_density(psi)
    return expect(psi, "ProjUp")
end

function dmrg_run_1d_chain(rb_over_a, delta_over_omega, config_dir="config", data_dir="data_dir"; outputlevel=0, write_when_maxdim_exceeds=5000, write_dir="SCRATCH")
    # Hamiltonian specification
    ham_config = JSON3.read(read(joinpath(config_dir, "ham_config.json")))
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
    dmrg_params = JSON3.read(read(joinpath(config_dir, "dmrg_config.json")))
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

function dmrg_run_square_lattice(rb_over_a, delta_over_omega, config_dir="config", data_dir="data_dir"; outputlevel=0, write_when_maxdim_exceeds=5000, write_dir="SCRATCH")
    # Hamiltonian specification
    ham_config = JSON3.read(read(joinpath(config_dir, "ham_config.json")))
    n_x = ham_config["n_x"]
    n_y = ham_config["n_y"]
    omega = 2π * ham_config["omega"]
    delta = delta_over_omega * omega

    Rb = (C6 / omega) ^ (1 / 6)
    a = Rb / rb_over_a

    hamiltonian_params = Dict("n_sites"=>n_sites, "omega"=>omega, "delta"=>delta, "a"=>a)
    
    if haskey(ham_config, "rr_cutoff")
        rr_cutoff = hamiltonian_params["rr_cutoff"]
        ham, sites = rydberg_ham_square_lattice(n_x, n_y, omega, delta, a; rr_cutoff=rr_cutoff)
    else
        ham, sites = rydberg_ham_square_lattice(n_x, n_y, omega, delta, a)
    end

    # DMRG specifications
    dmrg_params = JSON3.read(read(joinpath(config_dir, "dmrg_config.json")))
    nsweeps = dmrg_params["nsweeps"]
    maxdim = dmrg_params["maxdim"]
    noise = dmrg_params["noise"]
    cutoff = dmrg_params["cutoff"]
    krylovdim = dmrg_params["krylovdim"]

    # Run DMRG
    psi0 = randomMPS(sites, maxdim[1])
    energy_observer = DMRGObserver(energy_tol=dmrg_params["energy_tol"])

    write_path = joinpath(write_dir, "$(rb_over_a)_$(delta_over_omega)")
    mkpath(write_path)
    energy, psi = dmrg(ham, psi0; nsweeps, maxdim, cutoff, noise, eigsolve_krylovdim=krylovdim, observer=energy_observer, outputlevel=outputlevel, write_when_maxdim_exceeds=write_when_maxdim_exceeds, write_path=write_path)
    rm(write_path, recursive=true)

    # Save results
    params = Dict("hamiltonian_params" => hamiltonian_params, "dmrg_params" => dmrg_params)
    result = Dict("energy" => energy, "rydberg_density" => rydberg_density(psi))

    mkpath(data_dir)
    JSON3.write(joinpath(data_dir, "$(rb_over_a)-$(delta_over_omega)" * ".json"), Dict("result" => result, "params" => params))
end

end # module RydbergArrayDMRG
