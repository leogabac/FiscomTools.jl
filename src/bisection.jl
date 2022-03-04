

"""
    bisection(f, (a,b) ;TOL)

Approximates the root of `f(x)` in the interval from `a` and `b` using the bisection method.

# Examples
```julia-repl
julia> f(x) = x^2 - 1 ;
julia> bisection(f, (0,3), TOL = 1e-5 )
1.0000019073486328

```

"""
function bisection(f, ansatz::Tuple{Number,Number}; TOL::Float64=0.01)
    (lower,upper) = ansatz
    mid = (lower + upper)/2;
    ea = 1; #Initialize
    
    # We do the first iteration
    if f(lower)*f(mid) < 0
        upper = mid
    elseif f(lower)*f(mid) > 0
        lower = mid
    else
        return mid
    end

    #Now we start doing the iterations until it is done
    while ea > TOL
        last = mid #We store the last approximation
        mid = (lower + upper)/2
        ea = abs((mid-last)/mid);
        
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