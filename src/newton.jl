include("ndiff.jl")

# Simple Newton-Rhapson method that takes into account number of iterations
function newtonR(f,x::Number,step::Float64,n::Int64)::Float64
    for _ in 1:n
        x= x - f(x)/ndiff(f,x,step)
    end
    return x
end

function newtonR(f,x::Number,step::Int64,n::Int64)::Float64
    for _ in 1:n
        x= x - f(x)/ndiff(f,x,step)
    end
    return x
end

#Newton-Rhapson functions that take into account the error, i.e. whenever the fourth input is a Float64 type
function newtonR(f,x::Number,step::Float64,error::Float64)::Float64
    ea = 1; #initialize
    x = x= x - f(x)/ndiff(f,x,step) #first iteration, since it is a requirement to calculate the error
    while ea > error
        last = x
        x = x - f(x)/ndiff(f,x,step)
        ea = abs((x-last)/x)
    end
    return x
end

function newtonR(f,x::Number,step::Int64,error::Float64)::Float64
    ea = 1; #initialize
    x = x= x - f(x)/ndiff(f,x,step) #first iteration, since it is a requirement to calculate the error
    while ea > error
        last = x
        x = x - f(x)/ndiff(f,x,step)
        ea = abs((x-last)/x)
    end
    return x
end
