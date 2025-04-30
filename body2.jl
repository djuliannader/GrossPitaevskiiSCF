module body2
push!(LOAD_PATH, pwd())
using LinearAlgebra
export splitstep
import potential
import norm


function splitstep(l,n,b,k,epsilon1,epsilon2,maxit,ro)
   h = (2*l)/n
   dt = 10^(-2)
   tmax = maxit*dt 
   xvals = [i*h for i in (-n/2+1):(n/2-1)]
   dVop = [potential.V(xvals[i])  for i in 1:length(xvals)]
   Vop =  Array(Diagonal(dVop))
   psi0 = gaussian(xvals)
   psi = norm.normalizing(psi0,2*l/n)
   psiin = psi
   dk = 2*pi/(n*h)
   kvals1 = [i*dk for i in 0:(n/2-1)]
   kvals2 = [i*dk for i in (-n/2+1):-1]
   kvals=vcat(kvals1,kvals2) 
   kfac = Array(Diagonal(kvals))
   expK = exp(-0.5*kfac^2*dt)
   it=0
   dmu=1.0
   mu=1.0
   psi1 = psi0
   psi2 = psi0
   while (abs(dmu)>epsilon1) 
   #while (it<1000)
     dop1  = [abs2(psiin[i]) for i in 1:length(psiin)]
     psi1  = [exp(-0.5*(dVop[i]+b*dop1[i])*dt)*psiin[i] for i in 1:length(psiin)]
     psik  = fourierdis(psi1)
     psik  = [expK[i,i]*psik[i] for i in 1:length(psik)]
     psi2  = invfourierdis(psik)
     dop2  = [abs2(psi2[i]) for i in 1:length(psi2)]
     psi  = [exp(-0.5*(dVop[i]+b*dop2[i])*dt)*psi2[i] for i in 1:length(psi2)]
     psi   = norm.normalizing(psi,h)
     mu1   = mu
     mu    = (1/dt)*log(psiin[floor(Int,n/2)]/psi[floor(Int,n/2)])
     dmu   = abs(mu-mu1)/abs(mu)
     psiin = psi
     it=it+1
     #println("- iteration: ",it," Chem. Potential: ",dmu)
   end
   #psi = norm.normalizing(psi0,h)
   return [real(mu),it,psi,1]
end




function gaussian(list)
  n=length(list)
  gaussianlist = [exp(-list[i]^2) for i in 1:n]
  return gaussianlist
end

function fourierdis(list)
  n=length(list)
  fouriert=[0.0+0.0*im for i in 1:n]
  for s in 1:n
    for r in 1:n
    fouriert[s] = fouriert[s] + (1/n^(1/2))*list[r]*exp(2*pi*im*(r-1)*(s-1)/n)
    end
  end
  return fouriert
end

function invfourierdis(list)
  n=length(list)
  fouriert=[0.0+0.0*im for i in 1:n]
  for s in 1:n
    for r in 1:n
    fouriert[s] = fouriert[s] + (1/n^(1/2))*list[r]*exp(-2*pi*im*(r-1)*(s-1)/n)
    end
  end
  return fouriert
end


end