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
	 return [1,real(sumw),neg,fotoc,sumx,sump]
         end
         end


function wignerfpoint(psi,L,N,hbar,xinst,pinst)
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
	 sum=0.0+0*im
         ie=1
         ix = round(Int, (xinst + L)/d)
	 for k in imin:imax
	      y=x[k]
	      sum=sum+(exp(2*im*pinst*y/hbar))*conj.(psi[ix-iint+ie])*psi[ix+iint-ie+1]*d
	      ie=ie+1
	 end
	 w=sum/(hbar*pi)
	 return real(w)
         end

function wignerfnp(psi,L,N,hbar)
	 pi=acos(-1)
         imin=trunc(Int64,N/4)
	 iint=imin
         imax=length(psi)-imin
         Lp = L
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
	        sum=sum+(exp(2*im*pinst*y/hbar))*conj.(psi[i-iint+ie])*psi[i+iint-ie+1]*d
		ie=ie+1
	     end
	     w=sum/(hbar*pi)
	     sumw=sumw+  d*d*round(real(w),digits=16)
	     sumnw=sumnw + d*d*round(abs(real(w)),digits=16)
	     sumx2=sumx2 + d*d*round(real(w),digits=16)*(xinst*xinst)
             sumx= sumx  + d*d*round(real(w),digits=16)*(xinst)
             sump2= sump2 + d*d*round(real(w),digits=16)*(pinst*pinst)
             sump=sump+d*d*round(real(w),digits=16)*(pinst)
           end
	 end
         fotoc = (sump2-sump^2) + (sumx2-sumx^2)
         neg = real(sumnw) - real(sumw)
	 return [1,real(sumw),neg,fotoc,sumx,sump]
         end


function witnesses(psi,L,N,hbar,xav,pav)
	 pi=acos(-1)
         imin=trunc(Int64,N/4)
	 iint=imin
         imax=length(psi)-imin
         Lp = L
         d=2*L/N
	 x=[-L+i*d for i in 1:(N-1)] 
         eps=0.001
         sumw     = 0.0
         sumnw    = 0.0
         sumx2c   = 0.0
         sumx3c   = 0.0
         sumx4c   = 0.0
         sump2c   = 0.0
         sump3c   = 0.0
         sump4c   = 0.0
         sumpx    = 0.0
         sumppx   = 0.0
         sumpxx   = 0.0
         sumppxx  = 0.0
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
             sumw=sumw+  d*d*round(real(w),digits=16)
             sumnw=sumnw + d*d*round(abs(real(w)),digits=16)  
             sumx2c   = sumx2c  + d*d*round(real(w),digits=16)*(xinst-xav)^2
             sumx3c   = sumx3c  + d*d*round(real(w),digits=16)*(xinst-xav)^3       
             sumx4c   = sumx4c  + d*d*round(real(w),digits=16)*(xinst-xav)^4
             sump2c   = sump2c  + d*d*round(real(w),digits=16)*(pinst-pav)^2
             sump3c   = sump3c  + d*d*round(real(w),digits=16)*(pinst-pav)^3  
             sump4c   = sump4c  + d*d*round(real(w),digits=16)*(pinst-pav)^4
             sumpx    = sumpx   + d*d*round(real(w),digits=16)*((pinst-pav)*(xinst-xav))  
             sumppx   = sumppx  + d*d*round(real(w),digits=16)*((pinst-pav)^2*(xinst-xav))
             sumpxx   = sumpxx  + d*d*round(real(w),digits=16)*((pinst-pav)*(xinst-xav)^2)
             sumppxx  = sumppxx + d*d*round(real(w),digits=16)*((pinst-pav)^2*(xinst-xav)^2)      
           end
	 end
         xkurtosis = (sumx4c)/(3*sumx2c^2)
         pkurtosis = (sump4c)/(3*sump2c^2)
         QFIX = 4*(sumx2c)/(2*hbar)
         QFIP = 4*(sump2c)/(2*hbar)
         exp_n = (1/2)*(sump2c/hbar + pav^2 + (sumx2c/hbar + 1.0*xav^2)) - 1/2
         r1 = 4*pav^2*sump2c/hbar + 4*xav^2*sumx2c/hbar + 8*xav*pav*sumpx/hbar
         r2 = (sump4c - sump2c^2)/hbar^2 + (sumx4c - sumx2c^2)/hbar^2 + 2*(sumppxx - sumx2c*sump2c)/hbar^2
         r3 = 4*pav*(sump3c + sumpxx)/hbar^(3/2) + 4*xav*(sumppx + sumx3c)/hbar^(3/2)
         varn = (r1 + r2 + r3)/4 - 1/4
         neg = real(sumnw) - real(sumw)
         #println("test  ",real(sumnw)," ",neg)
	 return [real(sumnw),neg,xkurtosis,pkurtosis,QFIX,QFIP,4*varn/(4*exp_n),sumx2c + sump2c]
end

function witnesses2(psi,L,N,hbar)
	 pi=acos(-1)
         imin=trunc(Int64,N/4)
	 iint=imin
         imax=length(psi)-imin
         Lp = L
         d=2*L/N
	 x=[-L+i*d for i in 1:(N-1)] 
         eps=0.001
         sumw     = 0.0
         sumnw    = 0.0
         sumx2c   = 0.0
         sumx3c   = 0.0
         sumx4c   = 0.0
         sump2c   = 0.0
         sump3c   = 0.0
         sump4c   = 0.0
         sumpx    = 0.0
         sumppx   = 0.0
         sumpxx   = 0.0
         sumppxx  = 0.0
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
             sumw=sumw+  d*d*round(real(w),digits=16)
             sumnw=sumnw + d*d*round(abs(real(w)),digits=16)  
             sumx2c   = sumx2c  + d*d*round(real(w),digits=16)*(xinst-xav)^2
             sumx3c   = sumx3c  + d*d*round(real(w),digits=16)*(xinst-xav)^3       
             sumx4c   = sumx4c  + d*d*round(real(w),digits=16)*(xinst-xav)^4
             sump2c   = sump2c  + d*d*round(real(w),digits=16)*(pinst-pav)^2
             sump3c   = sump3c  + d*d*round(real(w),digits=16)*(pinst-pav)^3  
             sump4c   = sump4c  + d*d*round(real(w),digits=16)*(pinst-pav)^4
             sumpx    = sumpx   + d*d*round(real(w),digits=16)*((pinst-pav)*(xinst-xav))  
             sumppx   = sumppx  + d*d*round(real(w),digits=16)*((pinst-pav)^2*(xinst-xav))
             sumpxx   = sumpxx  + d*d*round(real(w),digits=16)*((pinst-pav)*(xinst-xav)^2)
             sumppxx  = sumppxx + d*d*round(real(w),digits=16)*((pinst-pav)^2*(xinst-xav)^2)      
           end
	 end
         xkurtosis = (sumx4c)/(3*sumx2c^2)
         pkurtosis = (sump4c)/(3*sump2c^2)
         QFIX = 4*(sumx2c)/(2*hbar)
         QFIP = 4*(sump2c)/(2*hbar)
         exp_n = (1/2)*(sump2c/hbar + pav^2 + (sumx2c/hbar + 1.0*xav^2)) - 1/2
         r1 = 4*pav^2*sump2c/hbar + 4*xav^2*sumx2c/hbar + 8*xav*pav*sumpx/hbar
         r2 = (sump4c - sump2c^2)/hbar^2 + (sumx4c - sumx2c^2)/hbar^2 + 2*(sumppxx - sumx2c*sump2c)/hbar^2
         r3 = 4*pav*(sump3c + sumpxx)/hbar^(3/2) + 4*xav*(sumppx + sumx3c)/hbar^(3/2)
         varn = (r1 + r2 + r3)/4 - 1/4
         neg = real(sumnw) - real(sumw)
         #println("test  ",real(sumnw)," ",neg)
	 return [real(sumnw),neg,xkurtosis,pkurtosis,QFIX,QFIP,4*varn/(4*exp_n),sumx2c + sump2c]
         end


end
