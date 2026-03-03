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


function witnesses(psi,L,N,hbar,xav,pav,theta)
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
         sumt1    = 0.0
         sumt2    = 0.0
         sumq1    = 0.0
         sumq2    = 0.0
         wmin     = 0.0
         xmin     = 0.0
         pmin     = hbar^(1/2)
         gwmin     = -1*10^(-4)
         gxmin     = 0.0
         gpmin     = hbar^(1/2)
         sumwnr   = 0.0
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
             sumt1    = sumt1 + d*d*round(real(w),digits=16)*(xinst*cos(theta) + pinst*sin(theta))^2
             sumt2    = sumt2 + d*d*round(real(w),digits=16)*(xinst*cos(theta) + pinst*sin(theta))
             sumq1    = sumq1 + d*d*round(real(w),digits=16)*(xinst*cos(theta-pi/2) + pinst*sin(theta-pi/2))^2
             sumq2    = sumq2 + d*d*round(real(w),digits=16)*(xinst*cos(theta-pi/2) + pinst*sin(theta-pi/2))
             if real(w) < gwmin
                 gwmin = real(w)
                 gxmin = xinst
                 gpmin = pinst
             end    
             if real(w) < 0
                   if (xinst^2 + pinst^2) < hbar
                       sumwnr = sumwnr + d*d*round(real(w),digits=16)
                       if real(w)<wmin
                           wmin = real(w)
                           xmin = xinst
                           pmin = pinst
                       end
                   end
               end    
           end
	 end
         xkurtosis = (sumx4c)/(3*sumx2c^2)
         pkurtosis = (sump4c)/(3*sump2c^2)
         QFIX = 4*(sumx2c)/(2*hbar)
         QFIP = 4*(sump2c)/(2*hbar)
         QFIT = 4*(sumt1 - sumt2^2)/(2*hbar)
         varq = (sumq1 - sumq2^2)
         exp_n = (1/2)*(sump2c/hbar + pav^2 + (sumx2c/hbar + 1.0*xav^2)) - 1/2
         r1 = 4*pav^2*sump2c/hbar + 4*xav^2*sumx2c/hbar + 8*xav*pav*sumpx/hbar
         r2 = (sump4c - sump2c^2)/hbar^2 + (sumx4c - sumx2c^2)/hbar^2 + 2*(sumppxx - sumx2c*sump2c)/hbar^2
         r3 = 4*pav*(sump3c + sumpxx)/hbar^(3/2) + 4*xav*(sumppx + sumx3c)/hbar^(3/2)
         varn = (r1 + r2 + r3)/4 - 1/4
         neg = real(sumnw) - real(sumw)
         #println("test  ",real(sumnw)," ",neg)
	 return [(1/2)*real(sumnw),neg,xkurtosis,pkurtosis,QFIX,QFIP,4*varn/(4*exp_n),sumx2c + sump2c,QFIT,hbar/(2*varq),gwmin,gxmin,gpmin,abs(sumwnr),wmin,xmin,pmin]
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

function Nozeroscut(psi,L,N,hbar,theta)
       int = 2*L/N
       Ni = trunc(Int, 1/int)
       xlist1 = [(1-int*i)*cos(theta+pi) for i in 0:(Ni-1)]
       xlist2 = [int*i*cos(theta) for i in 0:Ni]
       xlist=vcat(xlist1,xlist2)
       plist1 = [(1-int*i)*sin(theta+pi) for i in 0:(Ni-1)]
       plist2 = [int*i*sin(theta) for i in 0:Ni]
       plist=vcat(plist1,plist2)
       count = 0
       ep=10^(-3)
       for i in 1:(length(xlist)-1) 
           w1=wigner.wignerfpoint(psi,L,N,hbar,xlist[i],plist[i])
           w2=wigner.wignerfpoint(psi,L,N,hbar,xlist[i+1],plist[i+1])
       if (real(w1) > ep) && ( real(w2) < -ep)
         count = count + 1
       end
       if (real(w1) < -ep) && ( real(w2) > ep)
         count = count + 1
       end
       end
       return count
end

end
