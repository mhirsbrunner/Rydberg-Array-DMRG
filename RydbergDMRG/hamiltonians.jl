using ITensors

function build_ham(ham_config)
    lattice = lowercase(ham_config["lattice"])
    Rb = ham_config["Rb"]
    delta = ham_config["delta"]

    if cmp(lattice, "chain") * cmp(lattice, "1d_chain") * cmp(lattice, "1d") == 0
        n_sites = ham_config["n_sites"]

        if haskey(ham_config, "n_nns")
            ham, sites = ham_1d_chain(n_sites, Rb, delta, ham_config["n_nns"])
        else
            ham, sites = ham_1d_chain(n_sites, Rb, delta)
        end
    elseif cmp(lattice, "square") * cmp(lattice, "square_lattice") == 0
        n_x = ham_config["n_x"]
        n_y = ham_config["n_y"]

        if haskey(ham_config, "nn_cutoff")
            if haskey(ham_config, "y_periodic")
                ham, sites = ham_square_lattice(n_x, n_y, Rb, delta; y_periodic=ham_config["y_periodic"], nn_cutoff=ham_config["nn_cutoff"])
            else
                ham, sites = ham_square_lattice(n_x, n_y, Rb, delta; nn_cutoff=ham_config["nn_cutoff"])
            end
        elseif haskey(ham_config, "y_periodic")
            ham, sites = ham_square_lattice(n_x, n_y, Rb, delta; y_periodic=ham_config["y_periodic"])
        else
            ham, sites = ham_square_lattice(n_x, n_y, Rb, delta)
        end
    elseif cmp(lattice, "lieb") * cmp(lattice, "Lieb") == 0
        n_x = ham_config["n_x"]
        n_y = ham_config["n_y"]
        delta_local = ham_config["delta_local"]

        if haskey(ham_config, "nn_cutoff")
            if haskey(ham_config, "y_periodic")
                ham, sites = ham_lieb_lattice(n_x, n_y, Rb, delta, delta_local; y_periodic=ham_config["y_periodic"], nn_cutoff=ham_config["nn_cutoff"], bc_boundary=ham_config["bc_boundary"])
            else
                ham, sites = ham_lieb_lattice(n_x, n_y, Rb, delta, delta_local; nn_cutoff=ham_config["nn_cutoff"], bc_boundary=ham_config["bc_boundary"])
            end
        elseif haskey(ham_config, "y_periodic")
            ham, sites = ham_lieb_lattice(n_x, n_y, Rb, delta, delta_local; y_periodic=ham_config["y_periodic"], bc_boundary=ham_config["bc_boundary"])
        else
            ham, sites = ham_lieb_lattice(n_x, n_y, Rb, delta, delta_local, bc_boundary=ham_config["bc_boundary"])
        end
    end

    return ham, sites
end


############
# 1D Chain #
############

function ham_1d_chain(n_sites, Rb, delta, n_nns=4)
    if n_nns >= n_sites
        n_nns = n_sites - 1
    end

    sites = siteinds("S=1/2", n_sites)

    os = OpSum()

    for j=1:n_sites
        os += 1, "Sx", j
        os += -delta, "ProjUp", j
    end

    for r=1:n_nns
        for j=1:n_sites - r
            V = (Rb / r) ^ 6
            os += V, "ProjUp", j, "ProjUp", j + r
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

function ham_square_lattice(n_x, n_y, Rb, delta; yperiodic=true, nn_cutoff=2.1, pinning_field=0.0)
    n_sites = n_x * n_y
    
    lattice_sites = generate_sites_square_lattice(n_x, n_y)

    os = OpSum()

    # for ii = 1:2:n_y
    #     os += -pinning_field, "ProjUp", ii
    # end

    # os += -pinning_field, "ProjUp", 1

    for ii = 1:n_sites
        os += 1, "Sx", ii
        os += -delta, "ProjUp", ii
        os += -(1 + (-1)^ii) * pinning_field, "ProjUp", ii
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

            r = sqrt(x_disp ^ 2 + y_disp ^ 2)

            if r <= nn_cutoff
                V = (Rb / r) ^ 6
                os += V, "ProjUp", ind_1, "ProjUp", ind_2
            end
        end
    end

    sites = siteinds("S=1/2", n_sites)

    hamiltonian = MPO(os, sites)

    return hamiltonian, sites, os
end

# Currently unused
function site_to_index_square_lattice(site, n_x, n_y)
    if isodd(site[1])
        return n_y * (site[1] - 1) + site[2]
    else
        return n_y * (site[1] - 1) + 5 - site[2]
    end
end

##################
# Lieb Lattice #
##################

# function generate_sites_lieb_lattice(n_x, n_y, include_last_column=true)
#     # n_x and n_y are unit cells of the Lieb lattice, not rows and columns
#     sites = []

#     x_ind = 0
#     y_ind = 0

#     for ii in 1:n_x
#         x_ind += 1

#         for jj in 1:n_y
#             y_ind += 1
#             push!(sites, [x_ind, y_ind])
#             y_ind += 1
#             push!(sites, [x_ind, y_ind])
#         end

#         x_ind += 1
#         y_ind -= 1

#         for jj in 1:n_y
#             push!(sites, [x_ind, y_ind])
#             y_ind -= 2
#         end

#         y_ind += 1
#     end

#     if include_last_column
#         x_ind += 1

#         for jj in 1:n_y
#             y_ind += 1
#             push!(sites, [x_ind, y_ind])
#             y_ind += 1
#             push!(sites, [x_ind, y_ind])
#         end
#     end

#     return sites
# end

function generate_sites_lieb_lattice(n_x, n_y, y_periodic, bc_boundary=false)
    # n_x and n_y are unit cells of the Lieb lattice, not rows and columns
    atoms = []

    function add_a_site(x, y)
        push!(atoms, [2 * x, 2 * y])
    end

    function add_b_site(x, y)
        push!(atoms, [2 * x + 1, 2 * y])
    end

    function add_c_site(x, y)
        push!(atoms, [2 * x, 2 * y + 1])
    end
    

    for x in 0:n_x - 1
        for y in 0:n_y - 1
            add_a_site(x, y)
            add_b_site(x, y)
            add_c_site(x, y)

            if bc_boundary
                if x == 0
                    add_b_site(x - 1, y)
                end

                if y == 0 && !y_periodic
                    add_c_site(x, y - 1)
                end
            else
                if x == n_x - 1
                    add_a_site(x + 1, y)
                    add_c_site(x + 1, y)
                end

                if y == n_y - 1 && !y_periodic
                    add_a_site(x, y + 1)
                    add_b_site(x, y + 1)
                end
            end
        end
    end

    if !bc_boundary && !y_periodic
        add_a_site(n_x, n_y)
    end
             
    return atoms
end

function ham_lieb_lattice(n_x, n_y, Rb, delta, delta_local; y_periodic=true, nn_cutoff=2.1, bc_boundary=false)
    lattice_sites = generate_sites_lieb_lattice(n_x, n_y, y_periodic, bc_boundary)
    n_sites = length(lattice_sites)

    os = OpSum()

    for ii = 1:n_sites
        site = lattice_sites[ii]

        os += 1, "Sx", ii
        os += -delta, "ProjUp", ii
        
        if !(isodd(site[1]) && isodd(site[2]))
            os += delta_local, "ProjUp", ii
        end
    end

    for ind_1 = 1:n_sites
        for ind_2 = ind_1 + 1:n_sites
            site_1 = lattice_sites[ind_1]
            site_2 = lattice_sites[ind_2]

            x_disp = site_1[1] - site_2[1]
            temp_y_disp = site_1[2] - site_2[2]

            if y_periodic
                y_disp = min(abs(temp_y_disp), 2 * n_y - abs(temp_y_disp))
            else
                y_disp = temp_y_disp
            end

            r = sqrt(x_disp ^ 2 + y_disp ^ 2)

            if r <= nn_cutoff

                V = (Rb / r) ^ 6
                os += V, "ProjUp", ind_1, "ProjUp", ind_2
            end
        end
    end

    sites = siteinds("S=1/2", n_sites)

    hamiltonian = MPO(os, sites)

    return hamiltonian, sites, os
end