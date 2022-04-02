include("ndiff.jl")

"""
    newtonroot(f,x; step,TOL)

Computes the approximate root of `f(x)` with initial guess `x` using Newton's method.

The keyword (optional) arguments `step` defines the spacing for differentiation, and `TOL` the stopping criteria.

# Examples
```julia-repl
julia> f(x) = x^2 - 1;
julia> newtonroot(f,2)
1.0000000464611474

```
"""
function newtonroot(f,x::Real; step::Float64 = 0.01,TOL::Float64 = 0.01)
    ϵ = 1; #initialize
    x -=  f(x)/dv(f,x,h = step) #first iteration, since it is a requirement to calculate the error
    while ϵ > TOL
        last = x
        x -= f(x)/dv(f,x,h = step)
        ϵ = abs((x-last)/x)
    end
    return x
end