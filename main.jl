push!(LOAD_PATH, pwd())
import body
import norm
using Plots


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
 K7=readline(f)
 K8=readline(f)
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
wf=norm.normalizing(r[3],2L/N)
if K8=="True"
 println("u:",wf)
end
x=[-L+(2L/N)*i for  i in 1:(N-1)]
plot(x,wf,title="wave function")
xlabel!("x")
ylabel!("y(x)")
#if K8=="True"  
#   end
end



