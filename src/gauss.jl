
# Auxiliary functions
function swap_rows(A,i,j)
    A[i,:], A[j,:] = A[j,:],A[i,:] 
    return A
end

eye(n) = [ i==j ? 1.0 : 0.0 for i in 1:n, j in 1:n]


"""
    ref(A)

Computes the row echalon form of a square matrix ´A´. This is a subroutine of lgsolve(A).

For computation purposes, the matrix should be of type Matrix{Float64}.

# Examples
```julia-repl
julia> A = [1 3 3.0 4 5; 2 2 1 2 1; 3 2 3 3 7; 1 2 3 4 5; 7 8 9 4 5]
julia> ref(A)
5×5 Matrix{Float64}:
1.0  3.0  3.0   4.0         5.0
0.0  1.0  1.25  1.5         2.25
0.0  0.0  1.0   0.545455    2.81818
0.0  0.0  0.0   1.0        -1.55556
0.0  0.0  0.0   0.0       -23.3333

```
"""
function ref(A::Matrix{Float64})
    
    numRow, numCol = size(A); 
    
    for col in 1:numCol-1
        toone, tozero = eye(numRow), eye(numRow)
        # Is the pivot zero?
        if A[col,col] == 0.0 # Yes
            nonzero = findall(x->x!==0.0, A[:,col]) # Find all non-zero elements in the current column
            selected = nonzero[1] # We then select the first non-zero element in the column
            A = swap_rows(A,col,selected)
            toone[col,col] = 1/A[col,col]; # And make a new pivot
        else 
            toone[col,col] = 1/A[col,col]; 
        end
        tozero[col+1:end,col] = -A[col+1:end,col]
        
        A = tozero*toone*A

    end
    
    return A
end


"""
    lgsolve(A,b)

Computes the solution ´x´ of the matrix equation Ax = b.

# Examples
```julia-repl
julia> A = [1 3 3.0 4 5; 2 2 1 2 1; 3 2 3 3 7; 1 2 3 4 5; 7 8 9 4 5]
julia> b = [1,2,3,4,5.0];
julia> x = lgsolve(A,b)
julia> # From the LinearAlgebra package
julia> norm(A*x - b)
4.5288390936029406e-15

```
"""
function lgsolve(A::Matrix{Float64}, b::Vector{Float64})
    A = hcat(A,b);
    numRow, numCol = size(A); # Then retrieve some important variable names
    x = Vector{Float64}(undef,numRow)
    A = ref(A)
    
# Backwards solution
    x[end] = A[end,end];
    for k in reverse(1:length(x)-1)
        S=0; #initialize
        for m = 1:length(x)
            S += A[k,m]*x[m];
        end
        x[k] = A[k,end] - S; # isolate
    end
    return x
end