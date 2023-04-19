push!(LOAD_PATH, "/home/daniel/Documents/Bosons/GrossPitaevskii/Project")
import body


println("\r Gross-Pitaeskii 1D equation solver")
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
#------------------------------------

# printing information 
println("------------------------------------------------")
println("Space from -L to L with L=",L)
println("partitioning the space in N=",N, " subintervals")
println("Nonlinear term (beta)=",beta)

# calling function which performs selfconsistent method
@time begin
r=body.selfconsistent(L,N,beta)
end

# printing results
println("lambda = ",r[1])
iter= trunc(Int,r[2])
println("Convergence reached after ",iter, " iterations")
end


