using ITensors

C6 = 2π * .862690; # MHz μm^6

function build_ham(delta, a, ham_config)
    lattice = lowercase(ham_config["lattice"])
    omega = ham_config["omega"]

    if cmp(lattice, "chain") * cmp(lattice, "1d_chain") * cmp(lattice, "1d") == 0
        n_sites = ham_config["n_sites"]

        if haskey(ham_config, "n_nns")
            ham, sites = ham_1d_chain(n_sites, omega, delta, a, ham_config["n_nns"])
        else
            ham, sites = ham_1d_chain(n_sites, omega, delta, a)
        end
    elseif cmp(lattice, "square") * cmp(lattice, "square_lattice") == 0
        n_x = ham_config["n_x"]
        n_y = ham_config["n_y"]

        if haskey(ham_config, "nn_cutoff")
            if haskey(ham_config, "y_periodic")
                ham, sites = ham_square_lattice(n_x, n_y, omega, delta, a; y_periodic=ham_config["y_periodic"], nn_cutoff=ham_config["nn_cutoff"])
            else
                ham, sites = ham_square_lattice(n_x, n_y, omega, delta, a; nn_cutoff=ham_config["nn_cutoff"])
            end
        elseif haskey(ham_config, "y_periodic")
            ham, sites = ham_square_lattice(n_x, n_y, omega, delta, a; y_periodic=ham_config["y_periodic"])
        else
            ham, sites = ham_square_lattice(n_x, n_y, omega, delta, a)
        end
    end

    return ham, sites
end


############
# 1D Chain #
############

function ham_1d_chain(n_sites, omega, delta, a, n_nns=4)
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

##################
# Square Lattice #
##################

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

function ham_square_lattice(n_x, n_y, omega, delta, a; yperiodic=true, nn_cutoff=sqrt(2), pinning_field=0.000)
    n_sites = n_x * n_y
    
    lattice_sites = generate_sites_square_lattice(n_x, n_y)

    os = OpSum()

    for ii = 1:2:n_x
        os += -pinning_field, "ProjUp", ii
    end

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

            if r / a <= nn_cutoff
                os += V, "ProjUp", ind_1, "ProjUp", ind_2
            end
        end
    end

    sites = siteinds("S=1/2", n_sites)

    hamiltonian = MPO(os, sites)

    return hamiltonian, sites
end

# Currently unused
function site_to_index_square_lattice(site, n_x, n_y)
    if isodd(site[1])
        return n_y * (site[1] - 1) + site[2]
    else
        return n_y * (site[1] - 1) + 5 - site[2]
    end
end