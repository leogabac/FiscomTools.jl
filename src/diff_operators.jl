function d1hat(n::Int64, dx::Float64 )
    operator = zeros(n,n) # initialize
    for k in 1:n
        try operator[k, k+1] = 1. catch; Nothing end    
        try operator[k, k-1] = -1. catch; Nothing end
    end

    return operator./(2*dx)
end

function d2hat(n::Int64, dx::Float64 )
    operator = zeros(n,n)
    for k in 1:n
        operator[k,k] = -2.
        try operator[k, k+1] = 1. catch; Nothing end    
        try operator[k, k-1] = 1. catch; Nothing end
    end

    return operator./(dx^2)
end






