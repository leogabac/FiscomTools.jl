
"""
    trapz(f,a,b; h)

Computes the integral of `f(x)` from `a` to `b` using the trapezoidal rule.

The keyword (optional) argument defines the spacing between points. The default value is `h = 0.01`.

    trapz(x::Vector, y::Vector)

Computes the approximate integral using the discretized 1D mesh `x` and `y` using the trapezoidal rule. It is assumed that `x` is equally spaced.

# Examples
```julia-repl
julia> f(x) = cos(x) ;
julia> trapz(f,0,π/2, h = 1e-5 )
0.999999999939998

julia> x = collect(Float64,0:1e-5:π/2) ;
julia> y  = f.(x) ;
julia> trapz(x,y)
0.9999999999716525
```
"""
function trapz(f, a::Real, b::Real; h::Float64 = 0.01)
    Δ = (b-a)
    mid = sum( f(x) for x in (a+h:h:b-h) )
    I = Δ * ( f(a) + 2 * mid + f(b) ) / (2 * (Δ/h) );
    return I
end


# Integration of a table
function trapz(x::Vector, y::Vector)
    h, Δ = abs(x[2] - x[1]), x[end] - x[1]
    mid = sum( y[2:end-1] )
    I = Δ * (y[1] + 2 * mid + y[end]) / (2 * (Δ/h));
    return I
end