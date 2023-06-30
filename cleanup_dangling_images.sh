podman-hpc rmi -f $(podman-hpc images -f "dangling=true" -q)
