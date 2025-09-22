module main
push!(LOAD_PATH, pwd())
using Plots
include("modules/body.jl")
include("modules/body2.jl")
include("modules/norm.jl")
include("modules/energy.jl")
include("modules/wigner.jl")
include("modules/tunneling.jl")
include("modules/entropy.jl")
using .body
using .body2
using .norm
using .energy
using .wigner
using .tunneling
using .entropy






# Reading data from input file
#------------------------------------
open("input.dat") do f
 K1=readline(f)
 K2=readline(f)
 L = parse(Float64, K2)
 K3=readline(f)
 K4=readline(f)
 N = parse(Int64, K4)
 K5=readline(f)
 K6=readline(f)
 beta = parse(Float16, K6)
 K7=readline(f)
 K8=readline(f)
 k = parse(Int64, K8)
 K9=readline(f)
 K10=readline(f)
 K11=readline(f)
 K12=readline(f)
 ep1 = parse(Float32, K12)
 K13=readline(f)
 K14=readline(f)
 ep2 = parse(Float32, K14)
 K15=readline(f)
 K16=readline(f)
 mi = parse(Int64, K16)
 K17=readline(f)
 K18=readline(f)
 K19=readline(f)
 K20=readline(f)
 K21=readline(f)
 K22=readline(f)
 flagm = parse(Int64, K22)

 tt=r"([0-9])" 
 tpl = [parse(Int64,t.match) for t in eachmatch(tt, K20)]
 
#------------------------------------

# printing initial information
println("\r Gross-Pitaeskii 1D ")
if (flagm==2)
    println("\r Splitstep method")
else
    println("\r Selfconsistent method")
end
println("------------------------------------------------")
println("Space from -L to L with L=",L)
println("Partitioning the space in N=",N, " subintervals")
println("Nonlinear term (beta)=",beta)
println("------------------------------------------------")

# calling function which performs selfconsistent method
@time begin
if flagm==1
  r=body.selfconsistent(L,N,beta,k,ep1,ep2,mi)
end
if flagm==2
  r=body2.splitstep(L,N,beta,k,ep1,ep2,mi,-im)
end
end

# printing results
# calling function which normalize the wave function
wf=norm.normalizing(r[3],2L/N)
if r[4]==1
 iter= trunc(Int,r[2])
 println("-----------------------------------------------")
 # calling function which calculate the energy
 ener=energy.integratingenergy(wf,beta,L,N)
 println("*Convergence for the ",k," stationary state reached after ",iter," iterations*")
 entr=entropy.wehrlentropy(wf,L,N)
 println("Chemical potential = ",r[1])
 println("Energy             = ",real(ener))
 println("Wehrl Entropy      = ",real(entr))
 # calling function which calculate the derivative of the wave function at the center of coordinates
 der=norm.derivative2(wf,2L/N)
 println("Second derivative at the center of coordinates: ",der)
end

# Ploting the wave function
if r[4]==2
 println("Chemical potential = ",r[1])
 iter= trunc(Int,r[2])
 println("--->Failed to converge wave function after ",iter," iterations")
end

# calling function which calculate the classical turning points
 tpoints=tunneling.turnningpoints(wf,beta,L,N,r[1])
 println("Classical turning points: ",tpoints)

# calling function which transmision coefficient
 if K18=="True"
   if length(tpoints)>2
    tcoef=tunneling.wkbt(wf,tpoints[tpl[1]],tpoints[tpl[2]],beta,L,N,r[1])
    println("WKB transmission coeficient: ",tcoef)
   else
    println("*Make sure that there is an energy barrier at the energy of the stationary state")
   end
 end

# calling routine to calculate wigner function
   wig=wigner.wignerf(wf,L,N)


# printing the wave function
if K10=="True"
   d=2*L/N
   xx=[-L+i*d for i in 1:(N-1)]
   println("-------------")
   println("Go to file output/wavefunction.dat to see the wave function")
   println("-------------")
 open("output/wavefunction.dat","w") do io
 for i in 1:length(xx)
   println(io,xx[i]," ",real(r[3][i])," ",imag(r[3][i]))
 end
 end
end




end



end
