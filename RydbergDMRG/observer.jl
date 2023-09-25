using ITensors
mutable struct MyObserver <: AbstractObserver
    energy_tol::Float64
    entropy_tol::Float64
    trunc_tol::Float64
    min_sweeps::Int
    last_energy::Float64
    last_entropy::Float64
    last_m_s::Float64
    maxtruncerr::Float64
    MyObserver(energy_tol=0.0, entropy_tol=0.0, trunc_tol=0.0, min_sweeps=0) = new(energy_tol, entropy_tol, trunc_tol, min_sweeps, 1000.0, 1000.0, -1.0)
end

function ITensors.checkdone!(o::MyObserver;kwargs...)
    sweep = kwargs[:sweep]
    energy = kwargs[:energy]
    psi = kwargs[:psi]

    SvN = calculate_entranglement_entropy(psi)
    density = calculate_rydberg_density(psi)
    m_s = calculate_staggered_magnetization(psi)

    # println("SvN = $(SvN), m_s = $(m_s)")

    Ediff = abs(energy-o.last_energy)
    EEdiff = abs(SvN-o.last_entropy)
    m_sdiff = abs(m_s - o.last_m_s)

    if o.maxtruncerr <= o.trunc_tol && Ediff < o.energy_tol && EEdiff < o.entropy_tol && m_sdiff < 1e-3 && sweep > o.min_sweeps
        println("Stopping DMRG after sweep $sweep, energy error $(Ediff) < $(o.energy_tol), entropy error $(EEdiff) < $(o.entropy_tol),  and m_s error $(m_sdiff) <1e-3")
        return true
    end

    o.last_energy = energy
    o.last_entropy = SvN
    o.last_m_s = m_s
    o.maxtruncerr = -1
    
    flush(stdout)
    return false
end

function ITensors.measure!(o::MyObserver; kwargs...)
    energy = kwargs[:energy]
    sweep = kwargs[:sweep]
    bond = kwargs[:bond]
    outputlevel = kwargs[:outputlevel]
  
    spec = kwargs[:spec]
    o.maxtruncerr = max(o.maxtruncerr, spec.truncerr)

    if outputlevel > 2
      println("Sweep $sweep at bond $bond, the energy is $energy")
    end
end

function calculate_entranglement_entropy(psi)
    L = length(psi)
    b = div(L,2)

    psio = orthogonalize(psi, b)

    U,S,V = svd(psio[b], (linkind(psio, b-1), siteind(psio,b)))

    SvN = 0.0

    for n=1:dim(S, 1)
      p = S[n,n]^2
      SvN -= p * log(p)
    end

    return SvN
end