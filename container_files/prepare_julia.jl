using Pkg
Pkg.activate("RydbergArrayDMRG/")
Pkg.add.(["JSON3", "ITensors", "PackageCompiler", "ThreadPools", "Distributed", "Parallelism", "BenchmarkTools"])
Pkg.precompile()
Pkg.instantiate()

using PackageCompiler
create_sysimage(["JSON3", "ITensors", "RydbergArrayDMRG", "ThreadPools", "Distributed", "Parallelism", "BenchmarkTools"]; sysimage_path="RydbergArrayDMRGSysimage.so", precompile_execution_file="precompile_functions.jl")