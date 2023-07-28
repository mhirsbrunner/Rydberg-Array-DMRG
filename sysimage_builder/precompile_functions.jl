include("../RydbergDMRG/RydbergDMRG.jl")
using .RydbergDMRG

RydbergDMRG.run_phase_diagram_point(1.2, 1.2, "config_1d", "data";write_dir="data")
RydbergDMRG.run_phase_diagram_point(1.2, 1.3, "config_square", "data";write_dir="data")