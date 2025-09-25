module quenchdynamics
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("potential.jl")
include("norm.jl")
include("energy.jl")
using .potential
using .norm
using .energy
export dynamics        


function dynamics(l, n, b, th, nt, psi0, mphys; ħ=1.0)
    h = 2*l / n
    m = n - 1

    # Discrete Laplacian L (Dirichlet on ±l)
    d  = fill(-2/(h^2), m)
    du = fill( 1/(h^2), m-1)
    L  = Array(Tridiagonal(du, d, du))

    # Operators: kinetic T and diagonal potential A2
    T  = -(ħ^2/(2*mphys)) * L
    x  = -l .+ (1:m) .* h
    A2 = Array(Diagonal([potential.Vf(xj) for xj in x]))

    

    # (optional) higher-order combo coefficients you had
    c0, c1, c2 = 1.0, -1.0, 1.0

    psit = copy(psi0)
    psitt = copy(psi0)
    H = T + A2

    open("output/survivalamplitude.dat","w") do io
    for k in 1:Int(nt)
        println("- Time step ",k," of ",Int(nt))
        # Survival amplitude
        suramp = energy.integratingoverlap(psi0,psit,l,n)
        println(io,(k-1)*th," ",real(suramp)," ",imag(suramp))
        
        #diag(|psi|^2) at current state
        op1 = Array(Diagonal([abs2(psit[j]) for j in 1:m]))

        # predictor (Euler-type exponential step)
        psinp1til = exp((-1im/ħ)*th*(H + b*op1)) * psit
        psinp1til = norm.normalizing(psinp1til, h)

        # Midpoint estimate
        psinp12  = [(1/2)*(psit[i] + psinp1til[i]) for i in 1:m]
        psinp12  = norm.normalizing(psinp12, h)


        # Corrector (EM step with midpoint Hamiltonian)  
        op2  = Array(Diagonal([abs2(psinp12[j]) for j in 1:m]))
        psit = exp((-1im/ħ)*th*(H + b*op2)) * psit
        psit = norm.normalizing(psit, h)

            
        #psitt = exp((-1im/ħ)*th*H) * psitt
        #psitt = norm.normalizing(psitt, h)
    end
    end
    println("Go to file output/survivalamplitude.dat to see the Survival Amplitude")
    return psit
end


end
