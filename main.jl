push!(LOAD_PATH, pwd())
import body
import norm
import energy
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
println("partitioning the space in N=",N, " subintervals")
println("Nonlinear term (beta)=",beta)
println("------------------------------------------------")

# calling function which performs selfconsistent method
@time begin
r=body.selfconsistent(L,N,beta,k,ep1,ep2,mi)
end

# printing results
wf=norm.normalizing(r[3],2L/N)
x=[-L+(2L/N)*i for  i in 1:(N-1)]
if r[4]==1
 ener=energy.integratingenergy(wf,beta,L,N)
 println("chemical potential = ",r[1])
 println("energy             = ",ener)
 iter= trunc(Int,r[2])
 der=norm.derivative2(wf,2L/N)
 println("Second derivative at the origin of coordinates: ",der)
 if K10=="True"
   println("u:",wf)
 end
 println("Convergence for the ",k," state reached after ",iter," iterations")
end

if r[4]==2
 iter= trunc(Int,r[2])
 println("--->Failed to converge wave function after ",iter," iterations")
end
plot(x,wf,title="wave function")
 xlabel!("x")
 ylabel!("y(x)")
end



