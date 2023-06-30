FROM julia:latest

RUN apt-get update && \
    apt-get install -y gcc clang && \
    rm -rf /var/lib/apt/lists/*

COPY RydbergArrayDMRG ./RydbergArrayDMRG
COPY prepare_julia.jl .
COPY setup_config ./setup_config
COPY precompile_functions.jl .

RUN mkdir -p /shft
ENV JULIA_DEPOT_PATH=/shft

RUN julia prepare_julia.jl

COPY generate_phase_diagram.jl .