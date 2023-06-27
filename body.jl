module body
push!(LOAD_PATH, pwd())
using LinearAlgebra
export selfconsistent
import potential
import norm

function selfconsistent(l,n,b,k,epsilon1,epsilon2,maxit)
	 h=2*l/n
	 d=[-2/(h*h) for i in 1:(n-1)]
	 du=[1/(h*h) for i in 1:(n-2)]
	 A1=Array(Tridiagonal(du, d, du))
	 da2=[potential.V(-l+i*h) for i in 1:(n-1)]
	 A2=Array(Diagonal(da2))
	 Utemp=[]
	 Utemp=[1/((n-1)^(1/2)) for i in 1:(n-1)]
	 U=norm.normalizing(Utemp,2*l/n)
	 s=[0 for i in 1:(n-1)]
	 U=norm.normalizing(Utemp,2l/n)
	 A3=[]
	 da3=[]
	 H=[]
	 HV=[]
	 Ts=[]
	 lnew=2
	 it=0
	 s=[1 for i in 1:(n-1)]
	 while (sum(s)>epsilon2) && (it<maxit)
	 lold=10
	 while (abs(lold-lnew)>epsilon1) && (it<maxit)
	  it=it+1
	  da3=[U[i]*U[i] for i in 1:(n-1)]
	  A3=Array(Diagonal(da3))
	  H=(-1/2)*A1+A2+b*A3
	  Uold=[U[i] for i in 1:(n-1)]
	  HV=eigvecs(H)
	  Utemp=[HV[i,k] for i in 1:(n-1)]
	  #println("utemp:",Utemp)
	  #U=[HV[i,1] for i in 1:(n-1)]
	  U=norm.normalizing(Utemp,2*l/n)
	  lold=lnew
	  lnew=eigvals(H)[k]
	  s=[abs(U[i]-Uold[i]) for i in 1:(n-1)]
	 end
	 end
	 if it<maxit
	  return [lnew,it,U,1]
	 else 
	  return [lnew,it,U,2]
	 end
	 end 


end