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

# function delta_sweep(hamiltonian_params, delta_ax, dmrg_params, outputlevel=1)
#     # Unpack inputs
#     n_sites = hamiltonian_params["n_sites"]
#     omega = hamiltonian_params["omega"]
#     a = hamiltonian_params["a"]

#     nsweeps = dmrg_params["n_sweeps"]
#     maxdim = dmrg_params["max_dim"]
#     cutoff = dmrg_params["cutoff"]
#     krylov_dim = dmrg_params["krylov_dim"]
#     energy_tol = dmrg_params["energy_tol"]

#     energies = Array{Float32}(undef, length(delta_ax))
#     rydberg_densities = []

#     for (ii, delta) in enumerate(delta_ax)
#         if haskey(hamiltonian_params, "n_nns")
#             n_nns = hamiltonian_params["n_nns"]
#             ham, sites = rydberg_ham_1d(n_sites, omega, delta, a, n_nns)
#         else
#             ham, sites = rydberg_ham_1d(n_sites, omega, delta, a)
#         end

#         psi0 = randomMPS(sites, maxdim[1])
#         energy_observer = DMRGObserver(energy_tol=energy_tol)

#         energies[ii], psi = dmrg(ham, psi0; nsweeps, maxdim, cutoff, eigsolve_krylovdim=krylov_dim, observer=energy_observer, outputlevel=outputlevel)

#         push!(rydberg_densities, rydberg_density(psi))
#     end

#     return energies, rydberg_densities
# end
