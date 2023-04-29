push!(LOAD_PATH, pwd())
import body
import norm
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
 K9=readline(f)
 K10=readline(f)
 ep1 = parse(Float32, K10)
 K11=readline(f)
 K12=readline(f)
 ep2 = parse(Float32, K12)
#------------------------------------

# printing information 
println("------------------------------------------------")
println("Space from -L to L with L=",L)
println("partitioning the space in N=",N, " subintervals")
println("Nonlinear term (beta)=",beta)
println("------------------------------------------------")

# calling function which performs selfconsistent method
@time begin
r=body.selfconsistent(L,N,beta,ep1,ep2)
end

# printing results
wf=norm.normalizing(r[3],2L/N)
x=[-L+(2L/N)*i for  i in 1:(N-1)]
if r[4]==1
 println("lambda = ",r[1])
 iter= trunc(Int,r[2])
 if K8=="True"
   println("u:",wf)
 end
 println("Convergence reached after ",iter," iterations")
end
if r[4]==2
 iter= trunc(Int,r[2])
 println("--->Failed to converge wave function after ",iter," iterations")
end
plot(x,wf,title="wave function")
 xlabel!("x")
 ylabel!("y(x)")
end



