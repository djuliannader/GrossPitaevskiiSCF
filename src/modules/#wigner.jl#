module wigner
push!(LOAD_PATH, pwd())
include("body.jl")
include("norm.jl")
export wignerf


function wignerf(psi,L,N)
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
	 open("output/wignerfunction.dat","w") do io
	 sumw=0.0
	 sumnw=0.0
	 sumnx=0.0
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
	     println(io,xinst," ",pinst," ",round(real(w),digits=16))
	     sumw=sumw+d*d*w
	     sumnw=sumnw+d*d*abs(w)
	     sumnx=sumnx+d*d*w*(xinst*xinst)
           end
	 end
	 println("Go to file output/wignerfunction.dat to see the wigner function")
	 #println("Volume of the wigner function: ",real(sumw))
	 #println("Volume of the negative region: ",real(sumnw)-1)
	 #println("Expectation value <x^2> : ",real(sumnx))
	 return [1,real(sumw),real(sumnw)-1]
	 end
	 end


end
