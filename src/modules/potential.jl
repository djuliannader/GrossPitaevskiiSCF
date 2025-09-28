module potential
push!(LOAD_PATH, pwd())
export V
export Vf

# Initial potential
function V(str::AbstractString, q)
    str = replace(str, r"(\d)\s*([A-Za-z])" => s"\1*\2")
    # create a closure on the fly and apply it to q
    expr = Meta.parse("(x -> $str)($q)")
    return eval(expr)
end



end
