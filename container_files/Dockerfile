FROM julia:latest

# RUN apt-get update && \
#     apt-get install -y gcc clang && \
#     rm -rf /var/lib/apt/lists/*

RUN julia -e 'using Pkg; Pkg.add.(["JSON3", "ITensors", "ArgParse"]); Pkg.precompile(); Pkg.instantiate()'
# RUN julia -e 'using ITensors; ITensors.compile()'

COPY rydberg_arrays.jl .
COPY 1d_array_dmrg.jl .
COPY run_dmrg_point.sh .

RUN chmod +x run_dmrg_point.sh

ENTRYPOINT [ "/run_dmrg_point.sh" ]