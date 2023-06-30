@everywhere function func(a)
    RydbergArrayDMRG.dmrg_run_1d_chain(a...)
end

# BLAS.set_num_threads(2)
# ITensors.Strided.disable_threads()

phase_space_config = JSON3.read(read("config/phase_space_config.json"))

rb_over_a_ax = phase_space_config["rb_over_a_min"]:phase_space_config["rb_over_a_step"]:phase_space_config["rb_over_a_max"]

delta_over_omega_ax = phase_space_config["delta_over_omega_min"]:phase_space_config["delta_over_omega_step"]:phase_space_config["delta_over_omega_max"]

Parallelism.robust_pmap(func, Iterators.product(rb_over_a_ax, delta_over_omega_ax))

# @distributed for (rb_over_a, delta_over_omega) in Iterators.product(rb_over_a_ax, delta_over_omega_ax)
#     println("Starting Rb/a=$(rb_over_a), delta/omega=$(delta_over_omega)")
#     func(rb_over_a, delta_over_omega)
# end

# @sync for (rb_over_a, delta_over_omega) in Iterators.product(rb_over_a_ax, delta_over_omega_ax)
#     println("Starting Rb/a=$(rb_over_a), delta/omega=$(delta_over_omega)")
#     Threads.@spawn RydbergArrayDMRG.dmrg_run_1d_chain(rb_over_a, delta_over_omega)
# end

# qmap(RydbergArrayDMRG.dmrg_run_1d_chain, rb_over_a_ax, delta_over_omega_ax)
