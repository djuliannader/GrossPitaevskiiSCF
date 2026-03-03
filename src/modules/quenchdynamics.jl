module quenchdynamics
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("potential.jl")
include("norm.jl")
include("energy.jl")
include("wigner.jl")
include("body.jl")
using .potential
using .norm
using .body
using .wigner
using .energy
export dynamics        


function dynamics(l, n, b, th, nt, psi0, mphys, hbar, potf,theta)
    h = 2*l / n
    m = n - 1

    # Discrete Laplacian L (Dirichlet on ±l)
    d  = fill(-2/(h^2), m)
    du = fill( 1/(h^2), m-1)
    L  = Array(Tridiagonal(du, d, du))

    # Operators: kinetic T and diagonal potential A2
    T  = -(hbar^2/(2*mphys)) * L
    x  = -l .+ (1:m) .* h
    A2 = Array(Diagonal([potential.V(potf,xj) for xj in x]))

    

    # (optional) higher-order combo coefficients you had
    c0, c1, c2 = 1.0, -1.0, 1.0

    psit = copy(psi0)
    psitt = copy(psi0)
    H = T + A2

    open("output/survivalamplitude.dat","w") do io
    open("output/witnesses.dat","w") do io2      
    for k in 1:Int(nt)
        println("- Time step ",k," of ",Int(nt))
        # Survival amplitude
        suramp = energy.integratingoverlap(psi0,psit,l,n)
        wigs = wigner.wignerfnp(psit,l,n,hbar)
        wig = wigner.witnesses(psit,l,n,hbar,real(wigs[5]),real(wigs[6]),theta)
        nzero = wigner.Nozeroscut(psit,l,n,hbar,pi/2)
        println(io,(k-1)*th," ",real(suramp)," ",imag(suramp))
        println(io2,(k-1)*th," ",real(wig[1])," ",real(wig[2])," ",real(wig[5])," ",real(wig[6])," ",real(wig[7])," ",real(wig[9])," ",real(wig[10])," ",real(nzero[1])," ",real(wig[11])," ",real(wig[12])," ",real(wig[13])," ",real(wig[14])," ",real(wig[15])," ",real(wig[16])," ",real(wig[17]))
        
        #diag(|psi|^2) at current state
        op1 = Array(Diagonal([abs2(psit[j]) for j in 1:m]))

        # predictor (Euler-type exponential step)
        psinp1til = exp((-1im/hbar)*th*(H + b*op1)) * psit
        psinp1til = norm.normalizing(psinp1til, h)

        # Midpoint estimate
        psinp12  = [(1/2)*(psit[i] + psinp1til[i]) for i in 1:m]
        psinp12  = norm.normalizing(psinp12, h)


        # Corrector (EM step with midpoint Hamiltonian)  
        op2  = Array(Diagonal([abs2(psinp12[j]) for j in 1:m]))
        psit = exp((-1im/hbar)*th*(H + b*op2)) * psit
        psit = norm.normalizing(psit, h)

            
        #psitt = exp((-1im/ħ)*th*H) * psitt
        #psitt = norm.normalizing(psitt, h)
    end
    end
    end
    println("Go to file output/survivalamplitude.dat to see the Survival Amplitude")
    println("Go to file output/witnesses.dat to see some witnesses of the nonclassical dynamics")
    println("time | negativity  | X kurtosis | P kurtosis |  QFIX | QFIP | QFIN | FOTOC ")
    return psit
end


function complextimesp(wf0,L,N,beta,k,ep1,ep2,mi,mp,hbar,pot,tmax,gmax,tgint,Nlevels)
    enerlist = []
    clist = []
    for k in 1:Nlevels
        r=body.selfconsistentnp(L,N,beta,k,ep1,ep2,mi,mp,hbar,pot)
        wfn=norm.normalizing(r[3],2L/N)
        ev = energy.integratingenergy(wfn,beta,L,N,mp,hbar,pot)
        append!(enerlist,ev)
        cints = energy.integratingoverlap(wf0,wfn,L,N)
        append!(clist,abs2(cints))
    end
    open("output/survivalamplitudect.dat","w") do io
    for t in 0:tgint:tmax
    #println(t)        
        for g in -gmax:tgint:gmax
            zinst = 0.0 + 0.0*im
            for i in 1:length(enerlist)  
               zinst = zinst + clist[i]*exp(-im*(t+im*g)*enerlist[i]/hbar)
            end
        println(io,t," ",g," ",real(zinst)," ",imag(zinst))
            end
    end
    end
    return "done"
end
    
end
