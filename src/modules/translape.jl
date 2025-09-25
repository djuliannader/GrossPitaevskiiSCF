push!(LOAD_PATH, pwd())
include("body.jl")
include("norm.jl")
include("energy.jl")
include("tunneling,jl")
using .body
using .norm
using .energy
using .tunneling
using Plots


println("\r Gross-Pitaeskii 1D ")
println("\r Selfconsistent method")


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
#------------------------------------

# printing information 
println("------------------------------------------------")
println("Space from -L to L with L=",L)
println("Partitioning the space in N=",N, " subintervals")
println("Nonlinear term (beta)=",beta)
println("------------------------------------------------")

# calling function which performs selfconsistent method
@time begin
r1=body.selfconsistent(L,N,beta,k,ep1,ep2,mi)
kk=k+3
r2=body.selfconsistent(L,N,beta,kk,ep1,ep2,mi)
end

# printing results
wf1=norm.normalizing(r1[3],2L/N)
wf2=norm.normalizing(r2[3],2L/N)
x=[-L+(2L/N)*i for  i in 1:(N-1)]
if r1[4]==1
 ov=energy.integratingoverlap(wf1,wf2,L,N)
end

println("primer estado:",k)
println("segundo estado:",kk)
println("overlap:", ov)

end

