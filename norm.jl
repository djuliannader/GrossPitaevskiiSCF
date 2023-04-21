module norm
export normalizing

function normalizing(x::Vector{Float64},d)
    s=[i*i for i in x]
	 nfac=sum(d*s)
	 t=[i/(nfac)^(1/2) for i in x]
	 return t
	 end

end