# Richardson function
function ndiff(f,x::Number,h::Float64)
    D0 = (f(x+h)-f(x-h))/(2*h);
    D1 = ( f(x+(h/2)) - f(x-(h/2)) )/(h);
    TV = D1 + (D1-D0)/3
    return TV
end

#= 
# Richardson function
function ndiff(f,x::Number,h::Int64)
    D0 = (f(x+h)-f(x-h))/(2*h);
    D1 = ( f(x+(h/2)) - f(x-(h/2)) )/(h);
    TV = D1 + (D1-D0)/3
    return TV
end
=#

# Richardson discrete of type 1 (We get an exact value in x vector in which we want to differentiate)
function ndiff(x::Array{Float64,1}, y::Array{Float64,1}, val) #We get a x array, a y array and a val in x for which we want to get its derivative
    h = abs(x[2] - x[1]);
    ind = findall(p->p==val, x)[1]
    D0 = (y[ind+2]-y[ind-2])/(4*h); # The step in this case is double of the spacing
    D1 = (y[ind+1]-y[ind-1])/(2*h);
    TV = D1 + (D1-D0)/3
    return TV
end

# Richardson getting the whole bunch i.e. CDD of O(h^4) residues

#=
function ndiff(x::Array{Float64,1}, y::Array{Float64,1}) #We get a x array, a y array and a val in x for which we want to get its derivative
    h = abs(x[2] - x[1]);
    
    return TV
end
=#

# Forward func of linear order
function fdiff(f,x,h)
    d = (f(x+h)-f(x))/h;
    return d
end

# Forward discrete of type 1 (We get an exact value in x vector in which we want to differentiate)
function fdiff(x::Array{Float64,1}, y::Array{Float64,1},val)
    h = abs( x[2] - x[1] );
    ind = findall(p->p==val, x)[1]
    d = (y[ind+1] - y[ind])/h
    return d
end

# Forward discrete, with the whole bunch (MATLAB-like, with the whole vector)
function fdiff(x::Array{Float64,1}, y::Array{Float64,1})
    h = abs( x[2] - x[1] );
    d = (y[2:end] - y[1:end-1])/h
    return d
end

# Backwards func
function bdiff(f,x,h)
    d = (f(x)-f(x-h))/h;
    return d
end

# Backwards discrete of type 1 (We get an exact value in x vector in which we want to differentiate)
function bdiff(x::Array{Float64,1}, y::Array{Float64,1},val)
    h = abs( x[2] - x[1] );
    ind = findall(p->p==val, x)[1]
    d = (y[ind] - y[ind-1])/h
    return d
end

# Backwards discrete, with the whole bunch (MATLAB-like, with the whole vector)
function bdiff(x::Array{Float64,1}, y::Array{Float64,1})
    h = abs( x[2] - x[1] );
    d = (y[2:end] - y[1:end-1])/h
    return d
end
