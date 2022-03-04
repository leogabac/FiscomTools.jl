

# Auxiliary sub-routine
endpoint_prod(f,interval::Tuple) = f(interval[1])*f(interval[end])

"""
    incremental(f,min,max; N)

Computes all intervals containing roots between `min` and `max`. The return type is a vector of tuples.

The keyword (optional) argument defines the number of points in the discretization. The default value is `N = 5`.

This function works best along with the `bisection` function


# Examples
```julia-repl
julia> f(x) = x^2 - 5;
julia> incremental(f, -10,10)
2-element Vector{Tuple{Float64, Float64}}:
(-5.0, 0.0)
(0.0, 5.0)

```
"""
function incremental(f,min::Real,max::Real;N::Int64 = 5 )
    t = range(min,stop = max,length = N) # makes the domain
    intervals = [(t[i], t[i+1]) for i in 1:length(t)-1]
    possible = [ interval for interval in intervals if endpoint_prod(f,interval) < 0 ]
    return possible
end

"""
    bisection(f,a,b;TOL)

Approximates the root of `f(x)` in the interval between `a` and `b` using the bisection method.

The keyword (optional) argument `TOL` refers to the stop criteria. The default value is `TOL = 0.01`

# Examples
```julia-repl
julia> f(x) = x^2 - 1 ;
julia> bisection(f,0,3, TOL = 1e-5 )
1.0000019073486328

```

"""
function bisection(f, lower::Real ,upper::Real; TOL::Float64=0.01)
    mid = (lower + upper)/2;
    ε = 1; #Initialize
    
    # We do the first iteration
    if f(lower)*f(mid) < 0
        upper = mid
    elseif f(lower)*f(mid) > 0
        lower = mid
    else
        return mid
    end

    #Now we start doing the iterations until it is done
    while ε > TOL
        last = mid #We store the last approximation
        mid = (lower + upper)/2
        ε = abs((mid-last)/mid);
        
        if f(lower)*f(mid) < 0
            upper = mid
        elseif f(lower)*f(mid) > 0
            lower = mid
        else
            return mid
        end
        
    end

    return mid
end