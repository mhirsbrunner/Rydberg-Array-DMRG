using LinearAlgebra
@everywhere function func(a; outputlevel=0)
    ITensors.Strided.disable_threads()
    BLAS.set_num_threads(1)
    RydbergArrayDMRG.dmrg_run_1d_chain(a...; outputlevel=outputlevel)
end

phase_space_config = JSON3.read(read("config/phase_space_config.json"))

rb_over_a_ax = phase_space_config["rb_over_a_min"]:phase_space_config["rb_over_a_step"]:phase_space_config["rb_over_a_max"]

delta_over_omega_ax = phase_space_config["delta_over_omega_min"]:phase_space_config["delta_over_omega_step"]:phase_space_config["delta_over_omega_max"]

Parallelism.robust_pmap(func, Iterators.product(rb_over_a_ax, delta_over_omega_ax))