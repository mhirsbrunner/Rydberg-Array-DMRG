# using Pkg
# Pkg.add.(["JSON3", "ThreadPools", "Distributed", "Parallelism", "BenchmarkTools", "ThreadPinning", "PackageCompiler"])

using PackageCompiler

create_sysimage(["ParallelDataTransfer", "ThreadPinning", "ArgParse", "JSON3", "ThreadPools", "Distributed", "Parallelism", "BenchmarkTools"]; sysimage_path="sys_rydberg.so", precompile_execution_file="precompile_functions.jl", base_sysimage="sys_itensors.so")