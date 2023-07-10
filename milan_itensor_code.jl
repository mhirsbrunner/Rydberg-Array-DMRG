using ITensors
using Printf
using BSON
using JLD2
mutable struct MyObserver <: AbstractObserver
    energy_tol::Float64
    entropy_tol::Float64
    last_energy::Float64
    last_entropy::Float64
    MyObserver(energy_tol=0.0,entropy_tol=0.0) = new(energy_tol,entropy_tol,1000.0,1000.0)
end
function ITensors.checkdone!(o::MyObserver;kwargs...)
    sw = kwargs[:sweep]
    energy = kwargs[:energy]
    psi = kwargs[:psi]
    L = length(psi)
    b = div(L,2)
    psio = orthogonalize(psi, b)
    U,S,V = svd(psio[b], (linkind(psio, b-1), siteind(psio,b)))
    SvN = 0.0
    for n=1:dim(S, 1)
      p = S[n,n]^2
      SvN -= p * log(p)
    end
    println("SvN = $(SvN)")
    Ediff = abs(energy-o.last_energy)
    EEdiff = abs(SvN-o.last_entropy)
    if Ediff < o.energy_tol && EEdiff < o.entropy_tol && sw > 100
    println("Stopping DMRG after sweep $sw, energy error $(Ediff) < $(o.energy_tol) and entropy error $(EEdiff) < $(o.entropy_tol)")
    return true
    end
    # Otherwise, update last_energy and keep going
    o.last_energy = energy
    o.last_entropy = SvN
     
    if isfile(@sprintf "stop_dmrg_Energy%.3f" energy)
        rm(@sprintf "stop_dmrg_Energy%.3f" energy)
        return true
    end
    if isfile("stop_dmrg_all")
        rm("stop_dmrg_all")
        return true
    end
    
    return false
end
ITensors.space(::SiteType"HardCore") = 2
ITensors.state(::StateName"0", ::SiteType"HardCore") = [1.0, 0.0]
ITensors.state(::StateName"1", ::SiteType"HardCore") = [0.0, 1.0]
function ITensors.op!(Op::ITensor,
                        ::OpName"N",
                        ::SiteType"HardCore",
                        s::Index)
    Op[s'=>2,s=>2] = 1
end
function ITensors.op!(Op::ITensor,
                        ::OpName"Adag",
                        ::SiteType"HardCore",
                        s::Index)
    Op[s'=>1,s=>2] = 1
end
function ITensors.op!(Op::ITensor,
                        ::OpName"A",
                        ::SiteType"HardCore",
                        s::Index)
    Op[s'=>2,s=>1] = 1
end
        
function IndexConvert(a, b, nx, ny)
    if isodd(a)
        return Int((a-1)/2*3/2*ny + b)
    else 
        return Int((a-2)/2*3/2*ny+ ny+ ny/2+1 -(b+1)/2)
    end
end
function GenSites(nx, ny)
    atom_coordinate = [[1.0, 1.0]]
    for ii = 1:nx 
        if isodd(ii)
           for jj=1:ny
                push!(atom_coordinate, [ii, jj])
           end
        else
           for jj = 1: ny/2
                push!(atom_coordinate, [ii, 2*jj-1])
           end
        end  
    end
    popfirst!(atom_coordinate)
    return atom_coordinate 
end
function DMRG_Lieb(Nx, Ny, Om, Rb, Del, DelP, nsweeps, MaxD, trunc1)
    
    V0 = Rb^6*Om
    atom_coordinate = GenSites(Nx, Ny)
    Ntot = length(atom_coordinate)
    sites = siteinds("HardCore", Ntot)
    ampo = OpSum()   
    # onsite
    for j = 1:Ntot
        ampo += 1/2,"Adag",j
        ampo += 1/2,"A",j
    end
    for j = 1:Ntot
        ampo += -Del,"N",j
    end
    for ii=1: Ntot
        atom1 = atom_coordinate[ii]
        atom1x = atom1[1]
        atom1y = atom1[2]
        if (isodd(trunc(Int,atom1x))&&iseven(trunc(Int,atom1y))) || (iseven(trunc(Int,atom1x)))
            site_index1= IndexConvert(trunc(Int,atom1x), trunc(Int,atom1y), Nx, Ny)
            #@show ii, atom_coordinate[ii], site_index1
            ampo += -DelP,"N",site_index1
        end
    end
    for ii = 1: Ntot
        for jj = ii+1: Ntot
            
            atom1 = atom_coordinate[ii]
            atom1x = atom1[1]
            atom1y = atom1[2]
            site_index1= IndexConvert(atom1x, atom1y, Nx, Ny)
            atom2 = atom_coordinate[jj]
            atom2x = atom2[1]
            atom2y = atom2[2]
            site_index2= IndexConvert(atom2x, atom2y, Nx, Ny)
             
            distance12 = sqrt((atom1x-atom2x)^2 + (min(abs(atom1y-atom2y), Ny-abs(atom1y-atom2y)))^2)
            V12= V0/distance12^6
            if V12 >0.05
                ampo += V12,"N",site_index1,"N",site_index2
            end
        end
    end
 
    H = MPO(ampo, sites; cutoff=1E-15)
    println(maxdim.(H))
    psi0 = randomMPS(sites)
    sweeps = Sweeps(nsweeps)
    for j = 1:nsweeps
        stps0 = 50
        noise0 = 1E-12
        if j <= stps0
            sweeps.maxdim[j] = 10
            sweeps.noise[j] = 1E-1
        elseif 20*(div(j-stps0-1,10)+1) <= MaxD
            sweeps.maxdim[j] = 20*(div(j-stps0-1,10)+1)
            sweeps.noise[j] = 10.0^(-div(j-stps0-1,10)-2) >= noise0 ? 10.0^(-div(j-stps0-1,10)-2) : 0 #noise0
        else
            sweeps.maxdim[j] = MaxD
            sweeps.noise[j] = 0 #noise0
        end
    end
    setcutoff!(sweeps, trunc1)
    # Run the DMRG algorithm, returning energy 
    # (dominant eigenvalue) and optimized MPS
    
    obs = MyObserver(1E-9,1E-8)
    energy, psi = dmrg(H,psi0, sweeps, observer=obs, svd_alg = "qr_iteration", eigsolve_maxiter=1)
    println("Final energy = $energy")
#  
    b = div(Ntot,2)
    psio = orthogonalize(psi, b)
    Ut,St,Vt = svd(psio[b], (linkind(psio, b-1), siteind(psio,b)))
    SvN = 0.0
    RE = 0.0
    for n=1:dim(St, 1)
      p = St[n,n]^2
      SvN -= p * log(p)
      RE += p^2
    end
    RE = -log(RE) 
   
    
    NN_corr = correlation_matrix(psi, "N", "N")
    N_each= zeros(Ntot)
    for ii=1: Ntot
        
         N_each[ii]= NN_corr[ii, ii]
    end
    x1 = [1.0]
    y1 = [1.0]
    z1 = [0.0]
    
    for ii=1: Ntot
        atom1 = atom_coordinate[ii]
        atom1x = atom1[1]
        atom1y = atom1[2]
        site_index1= IndexConvert(atom1x, atom1y, Nx, Ny)
        push!(x1, atom1x)
        push!(y1, atom1y)
        push!(z1, N_each[site_index1])
    end
    popfirst!(x1)
    popfirst!(y1)
    popfirst!(z1)
   
    #fn= plot(scatter(x=x1,
    #y=y1, mode="markers",
    #marker=attr(size=26, color=z1, colorscale="Viridis", showscale=true)
    #))
    #savefig(fn, "Density_plot.png")
    return N_each, NN_corr, SvN, energy, RE
end
let
    # Parameters
    Nx = parse(Int64, ARGS[1])
    Ny = parse(Int64, ARGS[2])  # needs to be even for periodic boundary condition
    Rb = parse(Float64, ARGS[3]) # Rydberg blockade radius for V1 
    Di = parse(Float64, ARGS[4]) #detuning original
    Dn = parse(Int64, ARGS[5]) # number of detuning
    Dstep = parse(Float64, ARGS[6])
    Delta_vec = [Di+(i-1)*Dstep for i = 1:Dn]
    len_Delta = length(Delta_vec)
    Omega = 1.0
    # # this is the parameters for DMRG
    nsweeps = parse(Int64, ARGS[7]) #200
    MaxD = parse(Int64, ARGS[8]) #1000
    trunc1 = parse(Float64, ARGS[9])
    DelP=parse(Float64, ARGS[10])
    println("DelP=", round(DelP,digits=3))
    for bb= 1: len_Delta
        Del = Delta_vec[bb]
        path_str="//work//LAS//flint-lab//milan//Lieb//DelP//"
    if !isfile(path_str*"DMRG_Lieb_Lattice-NX$(Nx)_Ny$(Ny)_Rb_$(round(Rb,digits=3))_Del$(round(Del,digits=3))_DelP$(round(DelP,digits=3))_Omega$(Omega).jld2")
            N_each, NN_corr, SvN, energy, RE = DMRG_Lieb(Nx, Ny, Omega, Rb, Del, DelP, nsweeps, MaxD, trunc1)
            println(N_each)
            println(SvN)
            println("Step completed = $(bb)")
            jldsave(path_str*"DMRG_Lieb_Lattice-NX$(Nx)_Ny$(Ny)_Rb_$(round(Rb,digits=3))_Del$(round(Del,digits=3))_DelP$(round(DelP,digits=3))_Omega$(Omega).jld2"; Rb, Del, DelP, Omega, energy, SvN, RE, NN_corr, N_each)
            GC.gc()
    end
    end
end