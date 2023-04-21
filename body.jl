module body
push!(LOAD_PATH, pwd())
using LinearAlgebra
export selfconsistent
import potential

function selfconsistent(l,n,beta)
	 h=2*l/n
	 d=[-2/(h*h) for i in 1:(n-1)]
	 du=[1/(h*h) for i in 1:(n-2)]
	 A1=Array(Tridiagonal(du, d, du))
	 da2=[potential.V(-l+i*h) for i in 1:(n-1)]
	 A2=Array(Diagonal(da2))
	 U=[1/((n-1)^(1/2)) for i in 1:(n-1)]
	 A3=[]
	 da3=[]
	 H=[]
	 HV=[]
	 Ts=[]
	 lold=10
	 lnew=2
	 it=0
	 while (abs(lold-lnew)>1e-16)
	  it=it+1
	  da3=[U[i]^2 for i in 1:(n-1)]
	  A3=Array(Diagonal(da3))
	  H=(-1/2)*A1+A2+beta*A3
	  HV=eigvecs(H)
	  U=[HV[i,1] for i in 1:(n-1)]
	  lold=lnew
	  lnew=eigvals(H)[1]
	 end
	 return [lnew,it,[HV[i,1] for i in 1:(n-1)]]
	 end 


end