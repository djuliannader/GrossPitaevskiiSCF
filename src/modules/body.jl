module body
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("potential.jl")
include("norm.jl")
using .potential
using .norm
export selfconsistent        

function selfconsistent(l,n,b,k,epsilon1,epsilon2,maxit,mp,hbar,pot)
	 h=2*l/n
	 d=[-2/(h*h) for i in 1:(n-1)]
	 du=[1/(h*h) for i in 1:(n-2)]
	 A1=Array(Tridiagonal(du, d, du))
	 da2=[potential.V(pot,(-l+i*h)) for i in 1:(n-1)]
	 A2=Array(Diagonal(da2))
	 Utemp=[]
	 Utemp=[1/((n-1)^(1/2)) for i in 1:(n-1)]
	 U=norm.normalizing(Utemp,2*l/n)
	 #s=[0 for i in 1:(n-1)]
	 #U=norm.normalizing(Utemp,2*l/n)
	 A3=[]
	 da3=[]
	 H=[]
	 HV=[]
	 Ts=[]
	 lnew=2
	 it=0
	 s=[1 for i in 1:(n-1)]
	 while (abs(1-sum(s))>epsilon2) && (it<maxit)
	 lold=10
	 while (abs(lold-lnew)>epsilon1) && (it<maxit)
	  it=it+1
	  da3=[U[i]*U[i] for i in 1:(n-1)]
	  A3=Array(Diagonal(da3))
          #println("flag: ",mp)   
	  H=(-hbar^2/(2*mp))*A1+A2+b*A3
	  Uold=[U[i] for i in 1:(n-1)]
	  HV=eigvecs(H)
	  Utemp=[HV[i,k] for i in 1:(n-1)]
	  #println("utemp:",Utemp)
	  #U=[HV[i,1] for i in 1:(n-1)]
	  U=norm.normalizing(Utemp,2*l/n)
	  lold=lnew
	  lnew=eigvals(H)[k]
          #println(U)
	  s=[h*U[i]*Uold[i] for i in 1:(n-1)]
	  println("- iteration: ",it," Dev. Eigenvalue: ",abs(lnew-lold), " Dev. Eigenfunction ",abs(1-sum(s))," -")
	 end
	 end
	 if it<maxit
	  return [lnew,it,U,1]
	  #return [listeig,it,U,1]
	 else 
	  return [lnew,it,U,2]
	  #return [listeig,it,U,2]
	 end
	 end 


end
