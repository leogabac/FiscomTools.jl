#If an integer is passed down to the function, it means that the user is requesting an exact number of iterations.
function bisection(f,lower,upper, N::Int64)
    mid = 0
    for _ in 1:N
        
    mid = (lower + upper)/2
        
        if f(lower)*f(mid) < 0
            upper = mid
        else
            lower = mid
        end
        
    end
    return mid
end

#If then a Float is passed to the function, then it means it must mean an objective relative error.
function bisection(f,lower,upper, error::Float64)
    mid = (lower + upper)/2;
    ea = 1; #Initialize
#We do the first iteration, otherwise, there is no way to have a relative error to measure
    if f(lower)*f(mid) < 0
        upper = mid
    elseif f(lower)*f(mid) > 0
        lower = mid
    else
        return mid
    end
#Now we start doing the iterations until it is done
    while ea > error
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
