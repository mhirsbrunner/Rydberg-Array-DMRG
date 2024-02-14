using JSON3

function read_config(config_dir)
    ham_config = JSON3.read(read(joinpath(config_dir, "ham_config.json")), Dict)
    dmrg_config = JSON3.read(read(joinpath(config_dir, "dmrg_config.json")), Dict)

    return ham_config, dmrg_config
end

function store_results(data_dir, f_name, ham_config, dmrg_config, energy, ee, psi, rydberg_density)
    params = Dict("ham_config" => ham_config, "dmrg_config" => dmrg_config)
    
    results = Dict("energy" => energy, "entanglement_entropy" => ee, "rydberg_density" => rydberg_density)

    mkpath(data_dir)

    JSON3.write(joinpath(data_dir, f_name * ".json"), Dict("results" => results, "params" => params))
end

function store_results(data_dir, f_name, ham_config, dmrg_config, energy, ee, psi, rydberg_density, correlator)
    params = Dict("ham_config" => ham_config, "dmrg_config" => dmrg_config)
    
    results = Dict("energy" => energy, "entanglement_entropy" => ee, "rydberg_density" => rydberg_density, "correlator" => correlator)

    mkpath(data_dir)

    JSON3.write(joinpath(data_dir, f_name * ".json"), Dict("results" => results, "params" => params))
end

function store_results(data_dir, f_name, ham_config, dmrg_config, results)
    params = Dict("ham_config" => ham_config, "dmrg_config" => dmrg_config)
    
    mkpath(data_dir)

    JSON3.write(joinpath(data_dir, f_name * ".json"), Dict("results" => results, "params" => params))
end