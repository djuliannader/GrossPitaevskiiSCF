module potential
push!(LOAD_PATH, pwd())
export V

function V(x)
    V = -10*x^2+x^4
    #sig=1
    #L=10.0
    #V0=1.0
    #V =-V0*(exp(-(x-L/2)^2/sig^2)+exp(-(x+L/2)^2/sig^2))
    return V
    end

end
