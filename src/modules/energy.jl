module energy
push!(LOAD_PATH, pwd())
include("potential.jl")
using .potential
export integratingenergy
export integratingoverlap

function integratingenergy(u,b,LL,NN,mp)
    d=2*LL/NN
    x=[-LL+i*d for i in 1:(NN-1)]
    t1=[(1/(2*mp))*((u[i+1]-u[i])/d)^2*d for i in 1:(NN-2)]
    t2=[u[i]*u[i]*potential.V(x[i])*d for i in 1:(NN-2)]
    t3=[(b/2)*(u[i])^4*d for i in 1:(NN-2)]
    ener=sum(t1+t2+t3)
	 return ener
	 end

function integratingchempot(u,b,LL,NN)
    d=2*LL/NN
    x=[-LL+i*d for i in 1:(NN-1)]
    t1=[(1/2)*((u[i+1]-u[i])/d)^2*d for i in 1:(NN-2)]
    t2=[u[i]*u[i]*potential.V(x[i])*d for i in 1:(NN-2)]
    t3=[(b)*(u[i])^4*d for i in 1:(NN-2)]
    ener=sum(t1+t2+t3)
	 return ener
	 end

function integratingoverlap(u1,u2,LL,NN)
    d=2*LL/NN
    t2=[u1[i]*u2[i]*d for i in 1:(NN-1)]
    over=sum(t2)
	 return over
	 end


end
