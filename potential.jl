module potential
export V

function V(x)
    #V =-1.6874*x^2+x^4
    #V = -0.416452*x^2+0.0439106*x^4
    sig=1
    L=4*sig
    V =-(exp(-(x-L/2)^2/sig^2)+exp(-(x+L/2)^2/sig^2))
    return V
    end

end