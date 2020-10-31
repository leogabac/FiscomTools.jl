# Integration using a function
function trapz(f,a,b,h)
    x = collect(Float64,a:h:b);
    suma = 0
    for j ∈ 2:(length(x)-1)
        val = f(x[j])
        suma +=val
    end
    I = (b-a)*(f(x[1]) + 2*suma + f(x[end]))/(2*(length(x)-1));
    return I
end


# Integration of a table
function trapz(x::Array{Float64,1},y::Array{Float64,1})
    suma = 0;
    for j ∈ 2:(length(x)-1)
        suma +=y[j]
    end
    I = (x[end]-x[1])*(y[1] + 2*suma + y[end])/(2*(length(x)-1));
    return I
end