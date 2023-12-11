module RydbergDMRG

include("hamiltonians.jl")
include("observer.jl")
include("io.jl")

function calculate_rydberg_density(psi)
    return expect(psi, "ProjUp")
end

function staggered_magnetization_operator(psi)
    n_sites = length(psi)
    sites = siteinds(psi)

    M_N_os = OpSum()

    for ii=1:n_sites
        M_N_os += (-1) ^ ii / n_sites, "ProjUp", ii
        M_N_os += -1 / 2 * (-1) ^ ii / n_sites, "Id", ii
    end

    M_N = MPO(M_N_os, sites)

    return M_N
end

function colinear_order_parameter_operator(n_x, n_y, psi)
    sites = siteinds(psi)

    lattice_sites = generate_sites_lieb_lattice(n_x, n_y)
    n_sites = length(lattice_sites)

    os = OpSum()

    for ii=1:n_sites
        site = lattice_sites[ii]
        if isodd(site[1]) && !isodd(site[2])
            os += 1 / n_sites, "ProjUp", ii
        elseif !isodd(site[1]) && isodd(site[2])
            os -= 1 / n_sites, "ProjUp", ii
        end
    end

    colinear_mpo = MPO(os, sites)

    return colinear_mpo
end

function calculate_order_parameter(psi, op_mpo)
    return inner(psi', op_mpo, psi)
end

function calculate_order_parameter_moments(psi, op_mpo)
    op_2 = apply(op_mpo, op_mpo)
    op_4 = apply(op_2, op_2)

    exp_op_2 = inner(psi', op_2, psi)
    exp_op_4 = inner(psi', op_4, psi)

    return exp_op_2, exp_op_4
end

function calculate_order_parameter_binder_cumulant(exp_op_2, exp_op_4)
    return 3 / 2 - exp_op_4 / (exp_op_2 ^ 2) / 2
end

function calculate_order_parameter_susceptibility(n_sites, op, exp_op_2)
    return n_sites * (exp_op_2 - op ^ 2)
end

function calculate_staggered_magnetization(psi)
    n_sites = length(psi)
    sites = siteinds(psi)

    M_N_os = OpSum()

    for ii=1:n_sites
        M_N_os += (-1) ^ ii / n_sites, "ProjUp", ii
        M_N_os += -1 / 2 * (-1) ^ ii / n_sites, "Id", ii
    end

    M_N = MPO(M_N_os, sites)

    exp_M_N = inner(psi', M_N, psi)
    return abs(exp_M_N)
end

function calculate_stag_mag_moments(psi)
    n_sites = length(psi)
    sites = siteinds(psi)

    M_N_os = OpSum()

    for ii=1:n_sites
        M_N_os += (-1) ^ ii / n_sites, "ProjUp", ii
        M_N_os += -1 / 2 * (-1) ^ ii / n_sites, "Id", ii
    end

    M_N = MPO(M_N_os, sites)
    
    M_N_2 = apply(M_N, M_N)
    M_N_4 = apply(M_N_2, M_N_2)

    exp_M_N_2 = inner(psi', M_N_2, psi)
    exp_M_N_4 = inner(psi', M_N_4, psi)

    return exp_M_N_2, exp_M_N_4
end

function calculate_binder_cumulant(M_N_2, M_N_4)
    return 3 / 2 - M_N_4 / (M_N_2 ^ 2) / 2
end

function calculate_susceptibility(n_sites, m_s, M_N_2)
    return n_sites * (M_N_2 - m_s ^ 2)
end

function run_phase_diagram_point(Rb, delta, config_dir, data_dir; outputlevel=0, write_dir=nothing)
    f_name = "$(Rb)_$(delta)"

    ham_config, dmrg_config = read_config(config_dir)

    if haskey(ham_config, "delta_local")
        delta_local = ham_config["delta_local"]
        f_name = "$(Rb)_$(delta)_$(delta_local)"
    end

    ham_config["Rb"] = Rb
    ham_config["delta"] = delta

    ham, sites = build_ham(ham_config)

    # Run DMRG
    max_sweeps = dmrg_config["max_sweeps"]
    min_sweeps = dmrg_config["min_sweeps"]
    maxdim = dmrg_config["maxdim"]
    cutoff = dmrg_config["cutoff"]
    eigsolve_krylovdim = dmrg_config["krylovdim"]
    noise = dmrg_config["noise"]
    e_tol = dmrg_config["energy_tol"]
    ee_tol = dmrg_config["entanglement_entropy_tol"]
    trunc_tol = dmrg_config["trunc_tol"]

    psi0 = randomMPS(sites, maxdim[1])

    observer = RydbergDMRG.MyObserver(e_tol, ee_tol, trunc_tol, min_sweeps)

    if write_dir != nothing
        write_path = joinpath(write_dir, f_name)
        mkpath(write_path)
    end
    
    energy, psi = dmrg(ham, psi0; nsweeps=max_sweeps, maxdim, cutoff, eigsolve_krylovdim, noise, observer=observer, outputlevel)
    rydberg_density = calculate_rydberg_density(psi)
    ee = calculate_entranglement_entropy(psi)

    if write_dir != nothing
        rm(write_path, recursive=true)
    end

    store_results(data_dir, f_name, ham_config, dmrg_config, energy, ee, psi, rydberg_density)

    return psi
end

function run_local_detuning_phase_diagram_point(Rb, delta, delta_local, config_dir, data_dir; outputlevel=0, write_dir="/pscratch/sd/m/mhirsbru")
    f_name = "$(Rb)_$(delta)_$(delta_local)"

    ham_config, dmrg_config = read_config(config_dir)

    ham_config["Rb"] = Rb
    ham_config["delta"] = delta
    ham_config["delta_local"] = delta_local

    ham, sites = build_ham(ham_config)

    # Run DMRG
    max_sweeps = dmrg_config["max_sweeps"]
    min_sweeps = dmrg_config["min_sweeps"]
    maxdim = dmrg_config["maxdim"]
    cutoff = dmrg_config["cutoff"]
    eigsolve_krylovdim = dmrg_config["krylovdim"]
    noise = dmrg_config["noise"]
    e_tol = dmrg_config["energy_tol"]
    ee_tol = dmrg_config["entanglement_entropy_tol"]
    trunc_tol = dmrg_config["trunc_tol"]

    psi0 = randomMPS(sites, maxdim[1])

    observer = RydbergDMRG.MyObserver(e_tol, ee_tol, trunc_tol, min_sweeps)

    write_path = joinpath(write_dir, f_name)
    mkpath(write_path)
    
    energy, psi = dmrg(ham, psi0; nsweeps=max_sweeps, maxdim, cutoff, eigsolve_krylovdim, noise, observer=observer, outputlevel)
    
    rydberg_density = calculate_rydberg_density(psi)
    ee = calculate_entranglement_entropy(psi)

    rm(write_path, recursive=true)

    store_results(data_dir, f_name, ham_config, dmrg_config, energy, ee, psi, rydberg_density)

    return psi
end

function calculate_gap_magnetization_ent_entropy(Rb, delta, config_dir, data_dir; outputlevel=0, write_dir="/pscratch/sd/m/mhirsbru")
    f_name = "$(Rb)_$(delta)"

    ham_config, dmrg_config = read_config(config_dir)

    if haskey(ham_config, "delta_local")
        delta_local = ham_config["delta_local"]
        f_name = "$(Rb)_$(delta)_$(delta_local)"
    end

    ham_config["Rb"] = Rb
    ham_config["delta"] = delta

    ham, sites = build_ham(ham_config)

    # Run DMRG
    max_sweeps = dmrg_config["max_sweeps"]
    min_sweeps = dmrg_config["min_sweeps"]
    maxdim = dmrg_config["maxdim"]
    cutoff = dmrg_config["cutoff"]
    eigsolve_krylovdim = dmrg_config["krylovdim"]
    noise = dmrg_config["noise"]
    e_tol = dmrg_config["energy_tol"]
    ee_tol = dmrg_config["entanglement_entropy_tol"]
    trunc_tol = dmrg_config["trunc_tol"]

    psi0 = randomMPS(sites, maxdim[1])

    observer = RydbergDMRG.MyObserver(e_tol, ee_tol, trunc_tol, min_sweeps)

    write_path = joinpath(write_dir, f_name)
    mkpath(write_path)
    
    energy_ground, psi_ground = dmrg(ham, psi0; nsweeps=max_sweeps, maxdim, cutoff, eigsolve_krylovdim, noise, observer=observer, outputlevel)

    energy_excited, psi_excited = dmrg(ham, [psi_ground], psi0; nsweeps=max_sweeps, maxdim, cutoff, eigsolve_krylovdim, noise, observer=observer, outputlevel)
    
    rydberg_density = calculate_rydberg_density(psi_ground)
    stag_mag = calculate_staggered_magnetization(psi_ground)
    ee = calculate_entranglement_entropy(psi_ground)
    energy_gap = energy_excited - energy_ground

    rm(write_path, recursive=true)

    results = Dict("rydberg_density" => rydberg_density, "stag_mag" => stag_mag, "entanglement_entropy" => ee, "groundstate_energy" => energy_ground, "energy_gap" => energy_gap)


    store_results(data_dir, f_name, ham_config, dmrg_config, results)
end

function calculate_finite_size_scaling_data(n_y, Rb, delta, config_dir, data_dir; outputlevel=0, write_dir="/pscratch/sd/m/mhirsbru", op_name=nothing)
    f_name = "$(n_y)_$(Rb)_$(delta)"

    ham_config, dmrg_config = read_config(config_dir)

    if haskey(ham_config, "delta_local")
        delta_local = ham_config["delta_local"]
        f_name = "$(n_y)_$(Rb)_$(delta)_$(delta_local)"
    end
    
    n_x = 2 * (n_y - 1)
    ham_config["n_y"] = n_y
    ham_config["n_x"] = n_x

    ham_config["Rb"] = Rb
    ham_config["delta"] = delta

    ham, sites = build_ham(ham_config)

    # Run DMRG
    max_sweeps = dmrg_config["max_sweeps"]
    min_sweeps = dmrg_config["min_sweeps"]
    maxdim = dmrg_config["maxdim"]
    cutoff = dmrg_config["cutoff"]
    eigsolve_krylovdim = dmrg_config["krylovdim"]
    noise = dmrg_config["noise"]
    e_tol = dmrg_config["energy_tol"]
    ee_tol = dmrg_config["entanglement_entropy_tol"]
    trunc_tol = dmrg_config["trunc_tol"]

    psi0 = randomMPS(sites, maxdim[1])

    write_path = joinpath(write_dir, f_name)
    mkpath(write_path)
    
    observer = RydbergDMRG.MyObserver(e_tol, ee_tol, trunc_tol, min_sweeps)

    energy, psi = dmrg(ham, psi0; nsweeps=max_sweeps, maxdim, cutoff, eigsolve_krylovdim, noise, observer=observer, outputlevel, svd_alg="divide_and_conquer")
    
    rydberg_density = calculate_rydberg_density(psi)

    if op_name == "colinear"
        op_mpo = colinear_order_parameter_operator(n_x, n_y, psi)
    elseif op_name == "star"
        op_mpo = star_order_parameter_operator(psi)
    else
        op_mpo = staggered_magnetization_operator(psi)
    end

    # m_s = calculate_staggered_magnetization(psi)
    # M_N_2, M_N_4 = calculate_stag_mag_moments(psi)
    # U_4 = calculate_binder_cumulant(M_N_2, M_N_4)
    # X_s = calculate_susceptibility(n_x * n_y, m_s, M_N_2)

    op = calculate_order_parameter(psi, op_mpo)
    exp_op_2, exp_op_4 = calculate_order_parameter_moments(psi, op_mpo)

    rm(write_path, recursive=true)

    results = Dict("energy" => energy, "rydberg_density" => rydberg_density, "op" => op, "exp_op_2" => exp_op_2, "exp_op_4" => exp_op_4)

    # println("m_s = $(m_s), X_s = $(X_s), U_4 = $(U_4)")

    store_results(data_dir, f_name, ham_config, dmrg_config, results)
end

end