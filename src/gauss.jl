

function swap_rows(A,i,j)
    A[i,:], A[j,:] = A[j,:],A[i,:] 
    return A
end

eye(n) = [ i==j ? 1.0 : 0.0 for i in 1:n, j in 1:n]

function redGauss(A::Matrix{Float64})
    
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

function lgsolve(A::Matrix{Float64})
    numRow, numCol = size(A); # Then retrieve some important variable names

    x = Vector{Float64}(undef,numRow)
    A = redGauss(A)
    
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


# ===== TESTING ===== #
A = [1 3 3.0 4 5; 2 2 1 2 1; 3 2 3 3 7; 1 2 3 4 5; 7 8 9 4 5]
b = [1,2,3,4,5.0]
A\b

x = lgsolve(hcat(A,b))

A*x-b