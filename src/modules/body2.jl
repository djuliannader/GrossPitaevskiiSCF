module body2
push!(LOAD_PATH, pwd())
using LinearAlgebra
export splitstep
include("potential.jl")
include("norm.jl")
include("energy.jl")
using .potential
using .norm
using .energy



function splitstep(l,n,b,k,epsilon1,epsilon2,maxit,ro,mp)
   h = (2*l)/n
   dt = 10^(-2)
   tmax = maxit*dt 
   xvals = [i*h for i in (-n/2+1):(n/2-1)]
   dVop = [potential.V(xvals[i])  for i in 1:length(xvals)]
   psi0 = gaussian(xvals)
   psi = norm.normalizing(psi0,2*l/n)
   psiin = psi
   dk = 2*pi/(n*h)
   kvals1 = [i*dk for i in 0:(n/2-1)]
   kvals2 = [i*dk for i in (-n/2+1):-1]
   kvals=vcat(kvals1,kvals2) 
   kfac = Array(Diagonal(kvals))
   expK = exp(-(0.5/mp)*kfac^2*dt)
   it=0
   dmu=1.0
   mu=1.0
   #psi2 = psi
   while (abs(dmu)>epsilon1)
   #while (it<10000)
     dop1  = [abs2(psiin[i]) for i in 1:length(psiin)]
     psi1  = [exp(-(0.5)*(dVop[i]+b*dop1[i])*dt)*psiin[i] for i in 1:length(psiin)]
     psik  = fourierdis(psi1)
     psik  = [expK[i,i]*psik[i] for i in 1:length(psik)]
     psi2  = invfourierdis(psik)
     dop2  = [abs2(psi2[i]) for i in 1:length(psi2)]
     psi  = [exp(-(0.5)*(dVop[i]+b*dop2[i])*dt)*psi2[i] for i in 1:length(psi2)]
     psi   = norm.normalizing(psi,h)
     mu1   = mu
     mu    = (1/dt)*log(psiin[floor(Int,n/2)]/psi[floor(Int,n/2)])
     dmu   = abs(mu-mu1)/abs(mu)
     psiin = psi
     it=it+1
   end
   cp = energy.integratingchempot(psi,b,l,n)
   return [real(cp),it,psi,1]
end




function gaussian(list)
  n=length(list)
  gaussianlist = [exp(-list[i]^2) for i in 1:n]
  return gaussianlist
end


function fourierdis(list)
    n = length(list)
    fouriert = ComplexF64[0+0im for _ in 1:n]
    c = inv(sqrt(n))
    @inbounds for s in 0:n-1
        acc = 0.0 + 0.0im
        step = cis(2π*s/n)         # e^{i 2π s/n}, faster than exp(im*θ)
        ω = 1.0 + 0.0im            # current twiddle = step^r
        @inbounds @simd for r in 1:n
            acc += list[r] * ω
            ω *= step
        end
        fouriert[s+1] = c * acc
    end
    return fouriert
end

# function fourierdis(list)
#   n=length(list)
#   fouriert=[0.0+0.0*im for i in 1:n]
#   for s in 1:n
#     sumlist = [(1/n^(1/2))*list[r]*exp(2*pi*im*(r-1)*(s-1)/n) for r in 1:n]
#     fouriert[s] = sum(sumlist) 
#   end
#   return fouriert
# end

# function invfourierdis(list)
#   n=length(list)
#   fouriert=[0.0+0.0*im for i in 1:n]
#   for s in 1:n
#     sumlist = [(1/n^(1/2))*list[r]*exp(-2*pi*im*(r-1)*(s-1)/n) for r in 1:n]
#     fouriert[s] = sum(sumlist)
#   end
#   return fouriert
# end


function invfourierdis(list)
    n = length(list)
    fouriert = ComplexF64[0+0im for _ in 1:n]
    c = inv(sqrt(n))
    @inbounds for s in 0:n-1
        acc = 0.0 + 0.0im
        step = cis(-2π*s/n)      # e^{-i 2π s / n}
        ω = 1.0 + 0.0im
        @inbounds @simd for r in 1:n
            acc += list[r] * ω
            ω *= step            # reuse twiddle multiplicatively
        end
        fouriert[s+1] = c * acc
    end
    return fouriert
end


end
