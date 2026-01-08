module main
push!(LOAD_PATH, pwd())
#using Plots
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
 mp = parse(Float32, K10)
 K11=readline(f)
 K12=readline(f)
 hbar = parse(Float32, K12)   
 K13=readline(f)
 K14=readline(f)
 K15=readline(f)
 K16=readline(f)
 ep1 = parse(Float32, K16)
 K17=readline(f)
 K18=readline(f)
 ep2 = parse(Float32, K18)
 K19=readline(f)
 K20=readline(f)
 mi = parse(Int64, K20)
 K21=readline(f)
 K22=readline(f)
 K23=readline(f)
 K24=readline(f)
 K25=readline(f)
 K26=readline(f)
 flagm = parse(Int64, K26)
 K27=readline(f)
 K28=readline(f)
 pot = K28   

 tt=r"([0-9])" 
 tpl = [parse(Int64,t.match) for t in eachmatch(tt, K24)]
 
#------------------------------------

# printing initial information
println("\r Gross-Pitaeskii 1D ")
if (flagm==2)
    println("\r Splitstep method")
else
    println("\r Selfconsistent method")
end
println("------------------------------------------------")
println("Space from -L to L with L= ",L)
println("Partitioning the space in N= ",N, " subintervals")
println("Nonlinear term (beta)= ",beta)
println("mass (m)= ",mp)
println("hbar    = ",hbar)    
println("------------------------------------------------")

# calling function which performs selfconsistent method
@time begin
if flagm==1
  println("Selfconsistent field method selected")  
    r=body.selfconsistent(L,N,beta,k,ep1,ep2,mi,mp,hbar,pot)
    #r2=body.selfconsistent(L,N,beta,2,ep1,ep2,mi,mp,hbar,pot)
end
if flagm==2
  println("Split-time soliton method selected")  
  r=body2.splitstep(L,N,beta,k,ep1,ep2,mi,-im,mp,hbar,pot)
end
end

    println("here")

# printing results
# calling function which normalize the wave function
wf=norm.normalizing(r[3],2L/N)
#wf = wf=norm.normalizing((1/2^(1/2))*(r[3]+1.0*r2[3]),2L/N)   
if r[4]==1
 iter= trunc(Int,r[2])
 println("-----------------------------------------------")
 # calling function which calculate the energy
 ener=energy.integratingenergy(wf,beta,L,N,mp,hbar,pot)
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
 tpoints=tunneling.turnningpoints(wf,beta,L,N,r[1],pot)
 println("Classical turning points: ",tpoints)
 
# calling function which transmision coefficient
 if K18=="True"
   if length(tpoints)>2
    tcoef=tunneling.wkbt(wf,tpoints[tpl[1]],tpoints[tpl[2]],beta,L,N,r[1],pot)
    println("WKB transmission coeficient: ",tcoef)
   else
    println("*Make sure that there is an energy barrier at the energy of the stationary state")
   end
 end

# printing the wave function
if K14=="True"
   d=2*L/N
   xx=[(-L+i*d) for i in 1:(N-1)]
   println("-------------")
   println("Go to file output/wavefunction.dat to see the wave function")
   println("-------------")
 open("output/wavefunction.dat","w") do io
 for i in 1:length(xx)
   println(io,xx[i]," ",real(wf[i])," ",imag(wf[i]))
 end
 end
end

# calling routine to calculate wigner function
    wig=wigner.wignerf(wf,L,N,hbar)
    if wig[1]==1
       println("Volume of Wigner function of the state: ",wig[2])
       println("Negativity volume of the state :",wig[3])
       
    end






end



end
