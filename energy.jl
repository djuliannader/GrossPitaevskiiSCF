module energy
export integratingenergy
import potential

function integratingenergy(u::Vector{Float64},LL,NN)
    d=2*LL/NN
    x=[-LL+i*d for i in range(0,NN)]
	 return 1
	 end

end