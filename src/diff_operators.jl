
"""
    Lhat(n,h)

Computes the integral operator of size `n x n` with spacing `h`.

# Examples
```julia-repl
julia> Lhat(5,0.01);
5×5 Matrix{Float64}:
0.01  0.0   0.0   0.0   0.0
0.01  0.01  0.0   0.0   0.0
0.01  0.01  0.01  0.0   0.0
0.01  0.01  0.01  0.01  0.0
0.01  0.01  0.01  0.01  0.01
```
"""
function Lhat(n::Int64, h::Real)
    L = [ (row <= col) ? (1) : (0) for col in 1:n, row in 1:n]
    return h*L
end

"""
    D1hat(n,h)

Computes the (centered) first derivative operator of size `n x n` with spacing `h`.

# Examples
```julia-repl
julia> D1hat(5,0.01);
5×5 Matrix{Float64}:
  0.0   5.0   0.0   0.0  0.0
 -5.0   0.0   5.0   0.0  0.0
  0.0  -5.0   0.0   5.0  0.0
  0.0   0.0  -5.0   0.0  5.0
  0.0   0.0   0.0  -5.0  0.0
```
"""
function D1hat(n::Int64, h::Real )
    operator = zeros(n,n) # initialize
    for k in 1:n
        try operator[k, k+1] = 1. catch; Nothing end    
        try operator[k, k-1] = -1. catch; Nothing end
    end

    return operator./(2*h)
end

"""
    D2hat(n,h)

Computes the second derivative operator of size `n x n` with spacing `h`.

# Examples
```julia-repl
julia> D2hat(5,0.01);
5×5 Matrix{Float64}:
 -200.0   100.0     0.0     0.0     0.0
  100.0  -200.0   100.0     0.0     0.0
    0.0   100.0  -200.0   100.0     0.0
    0.0     0.0   100.0  -200.0   100.0
    0.0     0.0     0.0   100.0  -200.0
```
"""
function D2hat(n::Int64, h::Real )
    operator = zeros(n,n)
    for k in 1:n
        operator[k,k] = -2.
        try operator[k, k+1] = 1. catch; Nothing end    
        try operator[k, k-1] = 1. catch; Nothing end
    end

    return operator./(h^2)
end




