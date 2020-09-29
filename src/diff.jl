function diff(f,x::Float64,h::Float64)
    D0 = (f(x+h)-f(x-h))/(2*h);
    D1 = ( f(x+(h/2)) - f(x-(h/2)) )/(h);
    TV = D1 + (D1-D0)/3
    return TV
end

function diff(f,x::Int64,h::Float64)
    D0 = (f(x+h)-f(x-h))/(2*h);
    D1 = ( f(x+(h/2)) - f(x-(h/2)) )/(h);
    TV = D1 + (D1-D0)/3
    return TV
end

function diff(f,x::Float64,h::Int64)
    D0 = (f(x+h)-f(x-h))/(2*h);
    D1 = ( f(x+(h/2)) - f(x-(h/2)) )/(h);
    TV = D1 + (D1-D0)/3
    return TV
end

function diff(f,x::Int64,h::Int64)
    D0 = (f(x+h)-f(x-h))/(2*h);
    D1 = ( f(x+(h/2)) - f(x-(h/2)) )/(h);
    TV = D1 + (D1-D0)/3
    return TV
end