module main_quenchdynamics
push!(LOAD_PATH, pwd())
using Plots
include("modules/body.jl")
include("modules/body2.jl")
include("modules/norm.jl")
include("modules/energy.jl")
include("modules/wigner.jl")
include("modules/tunneling.jl")
include("modules/entropy.jl")
include("modules/quenchdynamics.jl")
using .body
using .body2
using .norm
using .energy
using .wigner
using .tunneling
using .entropy
using .quenchdynamics



# Reading data from input file
#------------------------------------
open("input_quenchdynamics.dat") do f
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
 ep1 = parse(Float32, K10)
 K11=readline(f)
 K12=readline(f)
 ep2 = parse(Float32, K12)
 K13=readline(f)
 K14=readline(f)
 mi = parse(Int64, K14)
 K15=readline(f)
 K16=readline(f)
 flagm = parse(Int64, K16)
 K17=readline(f)
 K18=readline(f)
 mp = parse(Float32, K18)  
 K19=readline(f)
 K20=readline(f)
 h = parse(Float32, K20)   
 K21=readline(f)
 K22=readline(f)
 nt = parse(Int64, K22)
 K23=readline(f)
 K24=readline(f)
 pot = K24
 K25=readline(f)
 K26=readline(f)
 potf = K26   

 
#------------------------------------

# printing initial information
println("\r Gross-Pitaeskii 1D ")
println("\r Quench dynamics")
println("------------------------------------------------")
println("Space from -L to L with L = ",L)
println("Partitioning the space in N = ",N, " subintervals")
println("Nonlinear term (beta) = ",beta)
println("Interval for each time step (h) = ",h)    
println("------------------------------------------------")

# 
@time begin

# ------------------ Obtaining the initial state ---------------------------    
if flagm==1
  r=body.selfconsistent(L,N,beta,k,ep1,ep2,mi,mp,pot)
end
if flagm==2
  r=body2.splitstep(L,N,beta,k,ep1,ep2,mi,-im,mp,pot)
end
if r[4]==1
    wf0=norm.normalizing(r[3],2L/N)
    println("----------------------------")
    println("|Initial state -> Obtained  |")
    println("---------------------------" )       
else
    println("----!!!...Failed to obtain initial state!!! try with the other numerical method")
    exit()
end
    eg = energy.integratingenergy(wf0,beta,L,N,mp,pot)
    println("Initial ground state energy: ",eg)
# -------------------------------------------------------------------------      

# ------------   Obtaining the dynamics and printing results ---------------------------      
wft = quenchdynamics.dynamics(L,N,beta,h,nt,wf0,mp,potf)
# calling routine to calculate wigner function
wig=wigner.wignerf(wft,L,N)
if wig[1]==1 
       println("Volume of Wigner function of the state: ",wig[2])
       println("Negativity volume of the state :",wig[3])
end    
d=2*L/N
xx=[-L+i*d for i in 1:(N-1)]
println("Go to file output/wavefunction.dat to see the wave function")
open("output/wavefunction.dat","w") do io
for i in 1:length(xx)
   println(io,xx[i]," ",real(wft[i])," ",imag(wft[i]))
end
end
# ------------   Obtaining the dynamics ----------------------------------      
    
end

end

end
