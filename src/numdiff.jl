
"""
    dv(f,x; h)

Computes the approximate derivative of `f(x)` evaluated at `x` using the Richardson scheme.

The keyword (optional) argument `h` defines the spacing. The default value is `h = 0.01`.

# Examples
```julia-repl
julia> f(x) = x^2 - 1;
julia> dv(f,2.5)
4.999999999999893

```
"""
function dv(f,x::Real; h::Float64 = 0.01)
    D0 = (f(x+h)-f(x-h))/(2*h);
    D1 = ( f(x+(h/2)) - f(x-(h/2)) )/(h);
    TV = D1 + (D1-D0)/3
    return TV
end

"""
    pdv(f,var,x...; h)

Computes the approximate partial derivative of `f(x1,...,xn)` with respecto to `x_var` evaluated at `x` using the Richardson scheme.

The keyword (optional) argument `h` defines the spacing. The default value is `h = 0.01`.

# Examples
```julia-repl
julia> f(x,y,z) = x + y^2 + 3*z^3;
julia> p = (0,1,1);
julia> pdv(f,3,p...)
8.999999999999904

```
"""
function pdv(f,var,x::Real ...;h = 0.01)
    htilde = Tuple(vcat([0 for _ in 1:var-1], h , [0 for _ in var+1:length(x)] ))
    htilde2 = htilde./2
    D0 = ( f( (x.+htilde)... ) - f( (x.-htilde)... ) )/(2*h);
    D1 = ( f( (x.+htilde2)... ) - f( (x.-htilde2)... ) )/(h);
    TV = D1 + (D1-D0)/3
    return TV    
end

"""
    gradient(f,x...; h)

Computes the approximate gradient of `f(x1,...,xn)` evaluated at `x`.

The keyword (optional) argument `h` defines the spacing. The default value is `h = 0.01`.

# Examples
```julia-repl
julia> f(x,y,z) = x + y^2 + 3*z^3;
julia> p = (0,1,1);
julia> gradient(f,p...)
3-element Vector{Float64}:
0.9999999999999787
1.9999999999999425
8.999999999999904

```
"""
function gradient(f,x::Real ...; h = 0.01)
    n = length(x)
    grad = [ pdv(f,var,x...;h) for var in 1:n ]
    return grad
end


# ===============================================
# Deprecated stuff
# ===============================================

# Richardson discrete of type 1 (We get an exact value in x vector in which we want to differentiate)
function numdiff(x::Array{Float64,1}, y::Array{Float64,1}, val) #We get a x array, a y array and a val in x for which we want to get its derivative
    h = abs(x[2] - x[1]);
    ind = findall(p->p==val, x)[1]
    D0 = (y[ind+2]-y[ind-2])/(4*h);
    D1 = (y[ind+1]-y[ind-1])/(2*h);
    TV = D1 + (D1-D0)/3
    return TV
end



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


