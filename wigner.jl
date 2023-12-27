module wigner
push!(LOAD_PATH, pwd())
export wignerf
import body
import norm


#L=20
#k=2
#N=400
#beta=0.0
#ep1=0.0000001
#ep2=0.001
#mi=100

#r=body.selfconsistent(L,N,beta,k,ep1,ep2,mi)
#wf=norm.normalizing(r[3],2L/N)



function wignerf(psi::Vector{Float64},L,N)
	 pi=acos(-1)
         imin=trunc(Int64,N/4)
	 iint=imin
	 imax=length(psi)-imin
	 d=2*L/N
	 x=[-L+i*d for i in 1:(N-1)]
	 eps=0.001
	 if (abs(psi[imin])+abs(psi[imax]))>eps
	   println("*Please, consider a larger domain to see the Wigner function*")
	   return "Done"
         end
	 open("wignerfunction.dat","w") do io
	 sumw=0.0
	 sumnw=0.0
         for i in imin:imax
	   xinst=x[i]
	   for j in imin:imax
	     pinst=x[j]
	     sum=0.0+0*im
	     ie=1
	     for k in imin:imax
	        y=x[k]
	        sum=sum+(exp(2*im*pinst*y))*psi[i-iint+ie]*psi[i+iint-ie+1]*d
		ie=ie+1
	     end
	     w=sum/(pi)
	     println(io,xinst," ",pinst," ",round(real(w),digits=16))
	     sumw=sumw+d*d*w
	     sumnw=sumnw+d*d*abs(w)
           end
	 end
	 println("Go to file wignerfunction.dat to see data for wigner function")
	 println("Volume of the wigner function:",real(sumw))
	 println("Vulume of the negative region:",real(sumnw)-1)
	 end
         return "Done"
	 end

#wig=wignerf(wf,L,N)

end