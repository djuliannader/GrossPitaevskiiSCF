module wigner
push!(LOAD_PATH, pwd())
include("norm.jl")
using .norm
export wignerf


function wignerf(psi,L,N,hbar)
	 pi=acos(-1)
         imin=trunc(Int64,N/4)
	 iint=imin
         imax=length(psi)-imin
         #Lp = L*hbar
         #d=2*Lp/N
         #x=[-Lp+i*d for i in 1:(N-1)]
         Lp = L
         d=2*L/N
	 x=[-L+i*d for i in 1:(N-1)] 
	 eps=0.001
	 if (abs(psi[imin])+abs(psi[imax]))>eps
	   println("*Please, consider a larger domain to see the Wigner function*")
	   return [2,"Done"]
         end
	 open("output/wignerfunction.dat","w") do io
	 sumw=0.0
	 sumnw=0.0
	 sumx2=0.0
         sumx=0.0
         sump2=0.0
         sump=0.0    
         for i in imin:imax
	   xinst=x[i]
	   for j in imin:imax
	     pinst=x[j]
	     sum=0.0+0*im
	     ie=1
	     for k in imin:imax
	        y=x[k]
	        sum=sum+(exp(2*im*pinst*y/hbar))*conj.(psi[i-iint+ie])*psi[i+iint-ie+1]*d
		ie=ie+1
	     end
	     w=sum/(hbar*pi)
	     println(io,xinst," ",pinst," ",round(real(w),digits=16))
	     sumw=sumw+  d*d*round(real(w),digits=16)
	     sumnw=sumnw + d*d*round(abs(real(w)),digits=16)
	     sumx2=sumx2 + d*d*round(real(w),digits=16)*(xinst*xinst)
             sumx= sumx  + d*d*round(real(w),digits=16)*(xinst)
             sump2= sump2 + d*d*round(real(w),digits=16)*(pinst*pinst)
             sump=sump+d*d*round(real(w),digits=16)*(pinst)
           end
	 end
	 println("Go to file output/wignerfunction.dat to see the wigner function")
         fotoc = (sump2-sump^2) + (sumx2-sumx^2)
         neg = real(sumnw) - real(sumw)
	 return [1,real(sumw),neg,fotoc]
	 end
         end

function wignerfnp(psi,L,N)
	 pi=acos(-1)
         imin=trunc(Int64,N/4)
	 iint=imin
	 imax=length(psi)-imin
	 d=2*L/N
	 x=[-L+i*d for i in 1:(N-1)]
	 eps=0.001
	 if (abs(psi[imin])+abs(psi[imax]))>eps
	   println("*Please, consider a larger domain to see the Wigner function*")
	   return [2,"Done"]
         end
	 sumw=0.0
	 sumnw=0.0
	 sumx2=0.0
         sumx=0.0
         sump2=0.0
         sump=0.0    
         for i in imin:imax
	   xinst=x[i]
	   for j in imin:imax
	     pinst=x[j]
	     sum=0.0+0*im
	     ie=1
	     for k in imin:imax
	        y=x[k]
	        sum=sum+(exp(2*im*pinst*y))*conj.(psi[i-iint+ie])*psi[i+iint-ie+1]*d
		ie=ie+1
	     end
	     w=sum/(pi)
	     #println(io,xinst," ",pinst," ",round(real(w),digits=16))
	     sumw=sumw+d*d*round(real(w),digits=16)
	     sumnw=sumnw + d*d*round(abs(real(w)),digits=16)
	     sumx2=sumx2 + d*d*round(real(w),digits=16)*(xinst*xinst)
             sumx= sumx  + d*d*round(real(w),digits=16)*(xinst)
             sump2= sump2 +d*d*round(real(w),digits=16)*(pinst*pinst)
             sump=sump+d*d*round(real(w),digits=16)*(pinst)
           end
	 end
	 #println("Go to file output/wignerfunction.dat to see the wigner function")
         fotoc = (sump2-sump^2) + (sumx2-sumx^2)
         neg = real(sumnw) - real(sumw)
	 return [1,real(sumw),neg,fotoc]
	 end


end
