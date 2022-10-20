
# Auxiliary functions
function simp13(f,a::Real,b::Real,h::Real)
    Δ = (b-a)
    N = floor(Δ/h)
    odd = sum( f(x) for x in a+h:2h:b-h ) 
    even = sum( f(x) for x in a+2h:2h:b-h ) 
    I = Δ * ( f(a) + 4*odd + 2*even + f(b) )/(3*N)
    return I
end

function simp38(f,a::Real,b::Real,h::Real)
    Δ = (b-a)
    return (Δ/8)*( f(a) + 3*f(a+h) + 3*f(a+2h) + f(b) )
end

"""
    intsimp(f,a,b,h)

Computes the integral of `f(x)` from `a` to `b` using the Simpson 1/3 and/or 3/8 rule. Assumes equal spacing ´h´

# Examples
```julia-repl
julia> f(x) = cos(x) ;
julia> intsimp(f,0,π/2,1e-5)
1.0000000000862006
```
"""
function intsimp(f,a::Real,b::Real,h::Real)
    N = floor((b-a)/h)
    if N % 2 == 0 #i.e. we have an even number of intervals
        return simp13(f,a,b,h)
    else #i.e. we have an odd number of intervals
        I1 = simp13(f,a,a + (N-3)*h ,h)
        I2 = simp38(f,a + (N-3)*h,b,h)
        return I1 + I2
    end
end


# Integration whenever a pair of vectors is given as a function, x and y
function intsimp(x::Array{Float64,1},y::Array{Float64,1})
    h = abs( x[2] - x[1] )

    if mod(length(x)-1,2) == 0 #i.e. we have an even number of intervals
        oddsum, evensum = 0,0 #initialize
        for  j in 2:2:(length(x) -1)
            oddsum += y[j]
        end
        for j in 3:2:(length(x) -1)
            evensum += y[j]
        end
        I = (x[end]-x[1])* ( y[1] + 4*oddsum + 2*evensum + y[end] )/(3*(length(x)-1))
        return I
    else #i.e. we have an odd number of intervals
        # First we need to calculate the simpson 1/3 up until... end-3
        # since we have N-3 points for simpson 1/3 and N-4 intervals
        oddsum, evensum = 0,0 #initialize
        for  j in 2:2:(length(x) - 4)
            oddsum += y[j]
        end
        for j in 3:2:(length(x) - 4)
            evensum += y[j]
        end
        I1 = (x[end-3]-x[1])* ( y[1] + 4*oddsum + 2*evensum + y[end-3] )/(3*(length(x)-4)) 

        I2 = (x[end]-x[end-3])/8*( y[end-3] + 3*y[end-2] + 3*y[end-1] + y[end] )

        return I1 + I2

    end
end


"""
    intsimp2(x,y,f)

Computes the 2D integral of `f(x,y)` in the rectangle delimited by `x` and `y` using `intsimp()`. The arguments must be given as a mesh.

# Examples
```julia-repl
julia> f(x,y) = sin(x)*cos(y) ;
julia> x = collect(-2:0.01:1) ;
julia> y = collect(-1:0.001:1)  ;
julia> feval = [f(xi,yi) for yi in y, xi in x] ;
julia> intsimp2(feval,x,y)
-1.609648403663146
```
"""
function intsimp2(f::Matrix,x::Vector,y::Vector)
    evaly = [intsimp(y,collect(col)) for col in eachcol(f)]
    return intsimp(x, evaly)
end
