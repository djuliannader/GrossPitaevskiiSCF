module entropy
push!(LOAD_PATH, pwd())
include("norm.jl")
include("body.jl")
using .norm
using .body
export wehrlentropy

function wehrlentropy(psi,LL,NN)
    d=2*LL/NN
    hbar=1.0
    x=[-LL+i*d for i in 1:(NN-1)]
    p=[-LL+i*d for i in 1:(NN-1)]
    vol=0.0
    for i in x
      for j in p
        alpha=coherentwf(i,j,hbar,LL,NN)
        tt=[conj(psi[k])*alpha[k]*d for k in 1:(NN-2)]
	qalpha=abs2(sum(tt))
	volal= qalpha*log(qalpha)*d^2/(2*pi)
	vol = vol - volal
      end
    end
    entropy=vol
    return entropy
end


function coherentwf(xav,pav,hbar,LL,NN)
         d=2*LL/NN
         x=[-LL+i*d for i in 1:(NN-1)]
	 psi=[wvc(xav,pav,hbar,i) for i in x]
	 psin=norm.normalizing(psi,d)
	 return psin
end


function wvc(xp,pp,hbar,x)
   fx=(1/pi^(1/4))*exp(im*x*pp/hbar)*exp(-(x-xp)^2/2)
   return fx
end

#xav=0.5
#pav=1.0
#N=1000
#L=20
#hbar=1.0
#beta=0.0
#k=1
#ep1=10^(-8)
#ep2=10^(-3)
#mi=100
#testwf=coherentwf(xav,pav,hbar,L,N)
#r=body.selfconsistent(L,N,beta,k,ep1,ep2,mi)
#testwf=norm.normalizing(r[3],2L/N)

#ent = wehrlentropy(testwf,L,N)

#println(ent)

end
