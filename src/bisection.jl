

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