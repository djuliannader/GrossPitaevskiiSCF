module norm
export normalizing, derivative2

function normalizing(x::Vector{Float64},d)
    s=[i*i for i in x]
	 nfac=sum(d*s)
	 t=[i/(nfac)^(1/2) for i in x]
	 return t
	 end



function derivative2(x::Vector{Float64},h)
	n=trunc(Int,length(x)/2)+1
	der=(x[n+1]-2*x[n]+x[n-1])/h^2
        return der
	end

end